# TOOL annotate-cpdb.py: "Hypergeometric test for ConsensusPathDB" (ConsensusPathDB created by Herwig et al. contains pathway information from over 20 publicly available databases, including Reactome, KEGG, BioCarta and HumanCyc. This service is provided by the Max Planck Institute for Molecular Genetics. NOTE! Supports any data with human, mouse or yeast gene symbols or UniProt identifiers.)
# INPUT input.tsv: input.tsv TYPE GENELIST 
# OUTPUT cpdb-pathways.html: cpdb-pathways.html 
# OUTPUT cpdb-pathways.tsv: cpdb-pathways.tsv 
# OUTPUT cpdb-genes.tsv: cpdb-genes.tsv 
# PARAMETER p_value_threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER input_type: input.type TYPE [hgnc-symbol: "gene symbol", uniprot: uniprot] DEFAULT hgnc-symbol (What kind of identifiers input data contains)
# PARAMETER species: species TYPE [human: human, mouse: mouse, yeast: yeast] DEFAULT human (Select the species of the data)

import sys
import csv

def main():
        
    global input_type, p_value_threshold, species
        
    input = 'input.tsv'
    # only for debugging, usually these values are defined in SADL
    try:
        input_type, p_value_threshold, species
    except NameError:
        input_type = 'hgnc-symbol'
        p_value_threshold = 0.05
        species = 'human'
        print('Parameters are undefined, using default values: ' + input_type + ', ' + str(p_value_threshold) + ', ' + species)
    
    cpdb_tool(input, input_type, p_value_threshold, species)

def cpdb_tool(input, input_type, p_value_threshold, species):
    """Chipster ConsensusPathDB tool.
    
    The input file must have gene symbols or uniprot ids and the 'type' parameter
    must be either 'gene' or 'uniprot' correspondingly. At first, the genes are mapped to 
    cpdb ids and then these ids are sent to over representation analysis. Analysis results
    are written to tsv and html files where there is a separate row for each pathway.
    In addition, the original input file (with genes) is extended with the associated
    pathways.  
    """
    
    if species == 'human':
        ws = 'http://cpdb.molgen.mpg.de/ws2/'
    elif species == 'mouse':
        ws = 'http://cpdb.molgen.mpg.de/ws2mouse/'
    elif species == 'yeast':
        ws = 'http://cpdb.molgen.mpg.de/ws2yeast/'    
    
    loc = cpdbLocator(ws)
    # Print soap messages to standard out
    # proxy = loc.getcpdb_portType(tracefile=sys.stdout)
    global proxy
    proxy = loc.getcpdb_portType()
    
    
    input_header = get_tsv_header(input)    
    genes = get_genes(input, input_header)
    
    
    id_to_gene_dict = query_gene_ids(genes, input_type)
            
    if len(id_to_gene_dict) == 0:
        raise Exception('None of the provided genes can be mapped to ' + \
                        'ConsensusPathDB identifiers. Check that input file has ' + \
                        'gene symbols or uniprot identifiers and the input type parameter is ' + \
                        'set correctly.')
    
    input_gene_ids = id_to_gene_dict.keys();
    
    def convert_ids_to_genes(ids, id_to_gene_dict):
        genes = []
        for id in ids:
            genes.append(id_to_gene_dict[id])
        return genes    
                                            
    pathways = query_over_representation_analysis(input_gene_ids, input_type, p_value_threshold)
    
    for pathway in pathways:                                                                                    
        pathway_genes = query_pathway_genes(pathway[PATHWAY_ID])
        
        if len(pathway_genes) != int (pathway[CORRECTED_GENE_COUNT]):
            raise Exception('Web service responds inconsistent number of pathway genes')
                          
        genes_mapped = convert_ids_to_genes(get_genes_mapped(pathway_genes, input_gene_ids), id_to_gene_dict)        
        
        if len(genes_mapped) != int (pathway[OVERLAPPING_GENE_COUNT]):
            raise Exception('Inconsistent number of overlapping genes')
    
        pathway[GENES_MAPPED] = genes_mapped    
    
    write_pathway_files(pathways, 'cpdb-pathways.tsv', 'cpdb-pathways.html')                
    write_extended_input_file(pathways, input, input_header, 'cpdb-genes.tsv')

def query_pathway_genes(pathway_id):
    """Get list of gene ids that are present in a pathway.
    
    Queries web service with the pathway id to get the contained genes.
    """    
    
    # result of getAvailableFsetTypesRequest():
    # 'P', manually curated pathways from pathway databases
    # 'N', interaction network neighborhood-based functional sets
    # 'G2', Gene Ontology-based sets, GO level 2
    # 'G3,' Gene Ontology-based sets, GO level 3
    # 'G4', Gene Ontology-based sets, GO level 4
    # 'G5', Gene Ontology-based sets, GO level 5
    # 'C', protein complex-based sets
    
    req = getCpdbIdsInFsetRequest()
    req._fsetId = pathway_id
    req._fsetType = 'P'
    req._entsetType = 'genes'
    response = proxy.getCpdbIdsInFset(req)
    result = response._cpdbIds
    return result

def get_gene(row, header):
    """Find a gene name on this row.
    
    Given a list of column values and a list of header names, this method returns the most suitable 
    value for gene name. If there is only one column, its value is returned. Otherwise columns 
    'symbol', ' ' and  'identifier' are searched in this order.  
    """
    
    gene_columns = ['symbol', ' ', 'identifier']
    
    if len(header) == 1:
        gene = row[0];
    else:
        for column in gene_columns:
            if column in header:
                gene = row[header.index(column)]
                break
    return gene

def create_links(pathway_name, db):
    """Create html links to pathway websites.            
    """    
    pathway_name = pathway_name.replace('_', ' ')
    encodedPathway = pathway_name.replace(' ', '+')    
    
    if db == 'Reactome':
        linked_name = '<a href=\"http://www.reactome.org/cgi-bin/search2?DB=gk_current&OPERATOR=ALL&QUERY=' + \
            encodedPathway + '&SPECIES=&SUBMIT=Go!\">' + pathway_name + '</a>'
            
    elif db == 'KEGG':
        pathwayWithoutOrganism = pathway_name[0:pathway_name.find('-')].strip();
        linked_name = '<a href=\"http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&serv=kegg&dbkey=kegg&keywords=' + \
            pathwayWithoutOrganism.replace(' ', '+') + '\">' + pathway_name + '</a>'

    elif db == 'PID':
        linked_name = '<a href=\"http://pid.nci.nih.gov/search/advanced_landing.shtml?' + \
            'what=graphic&svg=&jpg=true&xml=&biopax=&complex_uses=on&family_uses=on&degree=1&molecule=&' + \
            'pathway=' + encodedPathway + \
            '&macro_process=&source_id=5&evidence_code=NIL&evidence_code=IAE&evidence_code=IC&evidence_code=IDA&' + \
            'evidence_code=IFC&evidence_code=IGI&evidence_code=IMP&evidence_code=IOS&evidence_code=IPI&evidence_code=RCA&' + \
            'evidence_code=RGE&evidence_code=TAS&output-format=graphic&Submit=Go\">' + pathway_name + '</a>'

    elif db == 'HumanCyc':
        linked_name = '<a href=\"http://biocyc.org/HUMAN/substring-search?type=NIL&object=' + encodedPathway + '\">' + pathway_name + '</a>'
    
    else:
        linked_name = pathway_name

    return linked_name

def get_tsv_header(input):
    """Get a header of a tsv file.
    
    Column counts of the first (header) and second (content) line are compared. If the header has 
    one column less than the content, then a new column is added to the beginning of the returned list
    (convention in Chipster to handle R files).    
    """
    with open(input) as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)
        row = next(reader);     
     
    if len(header) + 1 == len(row):
        header = [' '] + header;
        
    return header

def get_genes(input, header):
    """Get a list of genes in a file.
    
    Reads the input file and creates a list of gene names or identifiers.
    """
    genes = []
       
    with open(input) as file:
        #skip header
        file.next()
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            gene = get_gene(row, header)                                 
            genes.append(gene)
            
    return genes

def query_gene_ids(genes, input_type):
    """Make a web service query to find cpdb gene identifiers for the given genes.
    
    Input is a list of gene symbols or uniprot identifiers. Results are returned
    as dict keys and the dict values contain the corresponding original genen symbol
    or uniprot identifier.
    """
    
    # result of getAvailableAccessionTypesRequest()
    # 'sgd-symbol', 'hgnc-id', 'sgd-id', 'entrez-gi', 'humancyc', 'intact', 'unigene', 'cygd', 'entrez-gene', 
    # 'hprd', 'mgi-symbol', 'mgi-id', 'hgnc-symbol', 'refseq', 'uniprot', 'ensembl', 'reactome', 'dip', 'pdb'
    
    # first map a list of uniprot IDs or gene symbols to cpdbIds
    req = mapAccessionNumbersRequest()
    req._accType = input_type # either hgnc-symbol or uniprot
    req._accNumbers = genes
    
    response = proxy.mapAccessionNumbers(req)
    result = zip(response._accNumber, response._cpdbId)
    
    id_to_gene_dict = dict()
    
    for r in result:
        if r[1]: # if mapping was found
            gene = r[0]
            id = r[1].split(',')[0]  ## here we take only one of the mappings for simplicity
            id_to_gene_dict[id] = gene

    return id_to_gene_dict

PATHWAY_NAME = 'pathway_name'
DB = 'db'
P_VALUE = 'p_value'
Q_VALUE = 'q_value'
GENE_COUNT = 'all_entities_num'
CORRECTED_GENE_COUNT = 'corrected_all_entities_num'
OVERLAPPING_GENE_COUNT = 'overlapping_entities_num'
PATHWAY_ID = 'fset_id'
URL = 'url'
# not in the original data, but will be added later
GENES_MAPPED = 'genes_mapped'

def query_over_representation_analysis(gene_ids, input_type, p_value_threshold):
    """Perform an over representation analysis for the given genes.
    
    Input is a list of cpdb gene ids. Returns a list dicts, i.e. a list of pathways.
    """
    req = overRepresentationAnalysisRequest()
    req._entityType = 'genes'
    req._fsetType = 'P'
    req._cpdbIdsFg = gene_ids
    ## req._cpdbIdsFg stays None, thus we have to set req._accType for the program to select the uniprot-specific background
    req._accType = input_type # either 'uniprot' or 'hgnc-symbol'
    # use user supplied parameter value
    req._pThreshold = p_value_threshold # default 0.05
    response = proxy.overRepresentationAnalysis(req)    
    
    cpdb_results = zip(response._name, response._details, response._overlappingEntitiesNum, response._allEntitiesNum, response._pValue, response._qValue)
    
    parsed_results = []
    
    for result in cpdb_results:
        pathway = dict()
        cpdb_name, cpdb_details, pathway[OVERLAPPING_GENE_COUNT], cpdb_all_entities_num, pathway[P_VALUE], pathway[Q_VALUE] = result
        pathway[PATHWAY_NAME], pathway[DB] = parse_name(cpdb_name)                
        pathway[GENE_COUNT], pathway[CORRECTED_GENE_COUNT] = parse_all_entities_num(cpdb_all_entities_num)
        pathway[PATHWAY_ID], pathway[URL] = parse_details(cpdb_details)
        parsed_results.append(pathway)        
                    
    return parsed_results

def parse_name(name):
    """Parse cpdb name field.
    """
    simple_name = name[0:name.rfind('(')]
    db = name[name.rfind('(') + 1:-1]
    
    return simple_name, db

def parse_all_entities_num(cpdb_all_entities_num):
    """Parse cpdb all entities num field.
    """
    all_entities_num = cpdb_all_entities_num[0:cpdb_all_entities_num.rfind('(')]
    corrected_all_entities_num = cpdb_all_entities_num[cpdb_all_entities_num.rfind('(') + 1:-1]
    
    return all_entities_num, corrected_all_entities_num

def parse_details(details):
    """Parse cpdb details field.
    """
    fsetIdString, sep, url = details.partition(';;')
    key, sep, fsetId = fsetIdString.partition(':')
    if not key == 'fsetId':
        raise Exeption("Can't find fsetId")
    return fsetId, url

def get_genes_mapped(pathway_genes, input_genes):
    """Returns an intersection of the two gene id lists. 
    """
    gene_ids_mapped = []
                     
    for cpdbId in pathway_genes:
        if cpdbId in input_genes:
            gene_ids_mapped.append(cpdbId)                    
    
    return gene_ids_mapped
    
def write_tsv(filename, header, output_rows):
    """Write a tab-separated file.
    """
    with open(filename, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(header)
        writer.writerows(output_rows)
        
def write_html(filename, header, output_rows):
    """Write the header and output_rows to file as html table.
    """
    
    html_header = '''
<html>
\t<html_header>
\t\t<style type="text/css">

h1 {
    padding: 8px;
    margin: 10px 0px 0px 24px;
    font-family: "Helvetica";
    font-weight: normal;
    font-size: 18;
}

table {
    margin: 10px 10px 0px 24px;
    border-width: 0px;
    border-spacing: 0px;
}

table th {
    padding: 8px;
    font-size: 11px;
    border-width: 0px 0px 0px 0px;
    border-bottom: 2px solid #0177b7;
    color: #0177b7;
    background: white;
    text-align: left;
}

table td {
    padding: 8px;
    font-size: 10px;
    border-width: 0px 0px 0x 0px;
    border-bottom: 1px solid gray;
    color: black;
}
\t\t</style>
\t\t<title>Over-representation analysis with ConsensusPathDB</title>
\t</html_header>
\t<body>
\t\t<h1>Over-representation analysis with ConsensusPathDB</h1>
'''        

    with open(filename, 'w') as file:

        file.write(html_header)
        file.write('\t\t<table>\n')
        # table header row
        file.write('\t\t\t<tr><th>')        
        file.write('</th><th>'.join(header))
        file.write('</th></tr>\n')
        # table content
        for row in output_rows:
            file.write('\t\t\t<tr><td>')        
            file.write('</td><td>'.join(row))
            file.write('</td></tr>\n')
            
        file.write('\t\t</table>\n')
        file.write('\t</body>\n')
        file.write('</html>\n')              
        
def write_pathway_files(pathways, tsv_filename, html_filename):
    """Write pathway information to tsv and html files.
    """
    # Filter columns and convert to str
    output_rows = []        
    for pathway in pathways:                                                                    
        genes_mapped_str = '; '.join(pathway[GENES_MAPPED])        
        output = (pathway[P_VALUE], pathway[OVERLAPPING_GENE_COUNT], pathway[CORRECTED_GENE_COUNT], pathway[PATHWAY_NAME], pathway[DB], genes_mapped_str)
        output_rows.append(output)

    pathway_header = ['p-value', 'Count', 'Size', 'Pathway', 'Database', 'Genes mapped']    

    write_tsv(tsv_filename, pathway_header, output_rows)

    # Create html links    
    html_rows = []
    for row in output_rows:
            html_row = list(row)
            html_row[3] = create_links(html_row[3], html_row[4])
            html_rows.append(html_row)
                                
    write_html(html_filename, pathway_header, html_rows)
    
def write_extended_input_file(pathways, input_filename, input_header, output_filename):
    """Copy the original input file and a new column for the pathways.
    """

    # a dict from gene name to list of pathway names 
    geneToPathwaysDict = dict()
        
    def add_pathway_to_genes(pathway_name, genes_mapped):
        """Store a pathway for all associated genes.
        """
        for gene in genes_mapped:
            if not gene in geneToPathwaysDict:
                geneToPathwaysDict[gene] = []
            geneToPathwaysDict[gene].append(pathway_name)        
            
    def get_pathways_of_gene(gene):
        """Get all pathways of the gene stored by the function add_pathway_to_genes().
        
        Pathways are converted to str using delimiter '; '.
        """
        if gene in geneToPathwaysDict:
            pathways_list = geneToPathwaysDict[gene]
            
            pathways_str = '; '.join(pathways_list)
            return pathways_str
        else:
            return ''
        
    # Input file has one gene on every row. Search pathways that mapped to this gene. 
    for pathway in pathways:                    
        add_pathway_to_genes(pathway[PATHWAY_NAME], pathway[GENES_MAPPED])
        
    # Read the input file, add a new pathway column and write the result to a file.
    with open(input_filename) as input:
        reader = csv.reader(input, delimiter='\t')
        with open(output_filename, 'w') as output:
            writer = csv.writer(output, delimiter='\t')
            #header
            row = next(reader)        
            row.append("pathway")
            writer.writerow(row)
            #content
            for row in reader:
                gene = get_gene(row, input_header)
                row.append(get_pathways_of_gene(gene))
                writer.writerow(row)

          
################################################## 
# cpdb_services_types.py 
# generated by ZSI.generate.wsdl2python
##################################################s

import ZSI
import ZSI.TCcompound
from ZSI.schema import LocalElementDeclaration, ElementDeclaration, TypeDefinition, GTD, GED

##############################
# targetNamespace
# cpdbns
##############################

class ns0:
    targetNamespace = "cpdbns"

    class getCpdbVersion_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getCpdbVersion"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getCpdbVersion_Dec.schema
            TClist = []
            kw["pname"] = ("cpdbns","getCpdbVersion")
            kw["aname"] = "_getCpdbVersion"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    return
            Holder.__name__ = "getCpdbVersion_Holder"
            self.pyclass = Holder

    class getCpdbVersionResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getCpdbVersionResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getCpdbVersionResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"cpdbVersion"), aname="_cpdbVersion", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getCpdbVersionResponse")
            kw["aname"] = "_getCpdbVersionResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._cpdbVersion = None
                    return
            Holder.__name__ = "getCpdbVersionResponse_Holder"
            self.pyclass = Holder

    class getAvailableAccessionTypes_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getAvailableAccessionTypes"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getAvailableAccessionTypes_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"entityType"), aname="_entityType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getAvailableAccessionTypes")
            kw["aname"] = "_getAvailableAccessionTypes"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._entityType = None
                    return
            Holder.__name__ = "getAvailableAccessionTypes_Holder"
            self.pyclass = Holder

    class getAvailableAccessionTypesResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getAvailableAccessionTypesResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getAvailableAccessionTypesResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"accType"), aname="_accType", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getAvailableAccessionTypesResponse")
            kw["aname"] = "_getAvailableAccessionTypesResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._accType = []
                    return
            Holder.__name__ = "getAvailableAccessionTypesResponse_Holder"
            self.pyclass = Holder

    class mapAccessionNumbers_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "mapAccessionNumbers"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.mapAccessionNumbers_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"accType"), aname="_accType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"accNumbers"), aname="_accNumbers", minOccurs=1, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","mapAccessionNumbers")
            kw["aname"] = "_mapAccessionNumbers"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._accType = None
                    self._accNumbers = []
                    return
            Holder.__name__ = "mapAccessionNumbers_Holder"
            self.pyclass = Holder

    class mapAccessionNumbersResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "mapAccessionNumbersResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.mapAccessionNumbersResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"accNumber"), aname="_accNumber", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"cpdbId"), aname="_cpdbId", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","mapAccessionNumbersResponse")
            kw["aname"] = "_mapAccessionNumbersResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._accNumber = []
                    self._cpdbId = []
                    return
            Holder.__name__ = "mapAccessionNumbersResponse_Holder"
            self.pyclass = Holder

    class getAvailableFsetTypes_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getAvailableFsetTypes"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getAvailableFsetTypes_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"entityType"), aname="_entityType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getAvailableFsetTypes")
            kw["aname"] = "_getAvailableFsetTypes"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._entityType = None
                    return
            Holder.__name__ = "getAvailableFsetTypes_Holder"
            self.pyclass = Holder

    class getAvailableFsetTypesResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getAvailableFsetTypesResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getAvailableFsetTypesResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"fsetType"), aname="_fsetType", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"description"), aname="_description", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getAvailableFsetTypesResponse")
            kw["aname"] = "_getAvailableFsetTypesResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._fsetType = []
                    self._description = []
                    return
            Holder.__name__ = "getAvailableFsetTypesResponse_Holder"
            self.pyclass = Holder

    class getDefaultBackgroundSize_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getDefaultBackgroundSize"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getDefaultBackgroundSize_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"fsetType"), aname="_fsetType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"accType"), aname="_accType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getDefaultBackgroundSize")
            kw["aname"] = "_getDefaultBackgroundSize"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._fsetType = None
                    self._accType = None
                    return
            Holder.__name__ = "getDefaultBackgroundSize_Holder"
            self.pyclass = Holder

    class getDefaultBackgroundSizeResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getDefaultBackgroundSizeResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getDefaultBackgroundSizeResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"bgSize"), aname="_bgSize", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getDefaultBackgroundSizeResponse")
            kw["aname"] = "_getDefaultBackgroundSizeResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._bgSize = None
                    return
            Holder.__name__ = "getDefaultBackgroundSizeResponse_Holder"
            self.pyclass = Holder

    class overRepresentationAnalysis_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "overRepresentationAnalysis"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.overRepresentationAnalysis_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"entityType"), aname="_entityType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"fsetType"), aname="_fsetType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"cpdbIdsFg"), aname="_cpdbIdsFg", minOccurs=1, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"cpdbIdsBg"), aname="_cpdbIdsBg", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"accType"), aname="_accType", minOccurs=0, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TCnumbers.FPfloat(pname=(ns,"pThreshold"), aname="_pThreshold", minOccurs=0, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","overRepresentationAnalysis")
            kw["aname"] = "_overRepresentationAnalysis"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._entityType = None
                    self._fsetType = None
                    self._cpdbIdsFg = []
                    self._cpdbIdsBg = []
                    self._accType = None
                    self._pThreshold = None
                    return
            Holder.__name__ = "overRepresentationAnalysis_Holder"
            self.pyclass = Holder

    class overRepresentationAnalysisResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "overRepresentationAnalysisResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.overRepresentationAnalysisResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"name"), aname="_name", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"details"), aname="_details", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"overlappingEntitiesNum"), aname="_overlappingEntitiesNum", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"allEntitiesNum"), aname="_allEntitiesNum", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"pValue"), aname="_pValue", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"qValue"), aname="_qValue", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","overRepresentationAnalysisResponse")
            kw["aname"] = "_overRepresentationAnalysisResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._name = []
                    self._details = []
                    self._overlappingEntitiesNum = []
                    self._allEntitiesNum = []
                    self._pValue = []
                    self._qValue = []
                    return
            Holder.__name__ = "overRepresentationAnalysisResponse_Holder"
            self.pyclass = Holder

    class enrichmentAnalysis_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "enrichmentAnalysis"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.enrichmentAnalysis_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"entityType"), aname="_entityType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"fsetType"), aname="_fsetType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"cpdbIdsMeasurements"), aname="_cpdbIdsMeasurements", minOccurs=1, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TCnumbers.FPfloat(pname=(ns,"pThreshold"), aname="_pThreshold", minOccurs=0, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","enrichmentAnalysis")
            kw["aname"] = "_enrichmentAnalysis"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._entityType = None
                    self._fsetType = None
                    self._cpdbIdsMeasurements = []
                    self._pThreshold = None
                    return
            Holder.__name__ = "enrichmentAnalysis_Holder"
            self.pyclass = Holder

    class enrichmentAnalysisResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "enrichmentAnalysisResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.enrichmentAnalysisResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"name"), aname="_name", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"details"), aname="_details", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"measuredEntitiesNum"), aname="_measuredEntitiesNum", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"allEntitiesNum"), aname="_allEntitiesNum", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"pValue"), aname="_pValue", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"qValue"), aname="_qValue", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","enrichmentAnalysisResponse")
            kw["aname"] = "_enrichmentAnalysisResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._name = []
                    self._details = []
                    self._measuredEntitiesNum = []
                    self._allEntitiesNum = []
                    self._pValue = []
                    self._qValue = []
                    return
            Holder.__name__ = "enrichmentAnalysisResponse_Holder"
            self.pyclass = Holder

    class getCpdbIdsInFset_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getCpdbIdsInFset"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getCpdbIdsInFset_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"fsetId"), aname="_fsetId", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"fsetType"), aname="_fsetType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded")), ZSI.TC.String(pname=(ns,"entsetType"), aname="_entsetType", minOccurs=1, maxOccurs=1, nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getCpdbIdsInFset")
            kw["aname"] = "_getCpdbIdsInFset"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._fsetId = None
                    self._fsetType = None
                    self._entsetType = None
                    return
            Holder.__name__ = "getCpdbIdsInFset_Holder"
            self.pyclass = Holder

    class getCpdbIdsInFsetResponse_Dec(ZSI.TCcompound.ComplexType, ElementDeclaration):
        literal = "getCpdbIdsInFsetResponse"
        schema = "cpdbns"
        def __init__(self, **kw):
            ns = ns0.getCpdbIdsInFsetResponse_Dec.schema
            TClist = [ZSI.TC.String(pname=(ns,"cpdbIds"), aname="_cpdbIds", minOccurs=0, maxOccurs="unbounded", nillable=False, typed=False, encoded=kw.get("encoded"))]
            kw["pname"] = ("cpdbns","getCpdbIdsInFsetResponse")
            kw["aname"] = "_getCpdbIdsInFsetResponse"
            self.attribute_typecode_dict = {}
            ZSI.TCcompound.ComplexType.__init__(self,None,TClist,inorder=0,**kw)
            class Holder:
                typecode = self
                def __init__(self):
                    # pyclass
                    self._cpdbIds = []
                    return
            Holder.__name__ = "getCpdbIdsInFsetResponse_Holder"
            self.pyclass = Holder

# end class ns0 (tns: cpdbns)

################################################## 
# cpdb_services.py 
# generated by ZSI.generate.wsdl2python
##################################################

import urlparse, types
from ZSI.TCcompound import ComplexType, Struct
from ZSI import client
import ZSI

# Locator
class cpdbLocator:
    def __init__(self, ws_address):
        self.cpdb_portType_address = ws_address
    def getcpdb_portTypeAddress(self):
        return self.cpdb_portType_address
    def getcpdb_portType(self, url=None, **kw):
        return cpdb_http_bindingSOAP(url or self.cpdb_portType_address, **kw)

# Methods
class cpdb_http_bindingSOAP:
    def __init__(self, url, **kw):
        kw.setdefault("readerclass", None)
        kw.setdefault("writerclass", None)
        # no resource properties
        self.binding = client.Binding(url=url, **kw)
        # no ws-addressing

    # op: getCpdbVersion
    def getCpdbVersion(self, request):
        if isinstance(request, getCpdbVersionRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#getCpdbVersion", **kw)
        # no output wsaction
        response = self.binding.Receive(getCpdbVersionResponse.typecode)
        return response

    # op: getAvailableAccessionTypes
    def getAvailableAccessionTypes(self, request):
        if isinstance(request, getAvailableAccessionTypesRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#getAvailableAccessionTypes", **kw)
        # no output wsaction
        response = self.binding.Receive(getAvailableAccessionTypesResponse.typecode)
        return response

    # op: mapAccessionNumbers
    def mapAccessionNumbers(self, request):
        if isinstance(request, mapAccessionNumbersRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#mapAccessionNumbers", **kw)
        # no output wsaction
        response = self.binding.Receive(mapAccessionNumbersResponse.typecode)
        return response

    # op: getAvailableFsetTypes
    def getAvailableFsetTypes(self, request):
        if isinstance(request, getAvailableFsetTypesRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#getAvailableFsetTypes", **kw)
        # no output wsaction
        response = self.binding.Receive(getAvailableFsetTypesResponse.typecode)
        return response

    # op: getDefaultBackgroundSize
    def getDefaultBackgroundSize(self, request):
        if isinstance(request, getDefaultBackgroundSizeRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#getDefaultBackgroundSize", **kw)
        # no output wsaction
        response = self.binding.Receive(getDefaultBackgroundSizeResponse.typecode)
        return response

    # op: overRepresentationAnalysis
    def overRepresentationAnalysis(self, request):
        if isinstance(request, overRepresentationAnalysisRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#overRepresentationAnalysis", **kw)
        # no output wsaction
        response = self.binding.Receive(overRepresentationAnalysisResponse.typecode)
        return response

    # op: enrichmentAnalysis
    def enrichmentAnalysis(self, request):
        if isinstance(request, enrichmentAnalysisRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#enrichmentAnalysis", **kw)
        # no output wsaction
        response = self.binding.Receive(enrichmentAnalysisResponse.typecode)
        return response

    # op: getCpdbIdsInFset
    def getCpdbIdsInFset(self, request):
        if isinstance(request, getCpdbIdsInFsetRequest) is False:
            raise TypeError("%s incorrect request type" % (request.__class__))
        kw = {}
        # no input wsaction
        self.binding.Send(None, None, request, soapaction="cpdbns#getCpdbIdsInFset", **kw)
        # no output wsaction
        response = self.binding.Receive(getCpdbIdsInFsetResponse.typecode)
        return response

getCpdbVersionRequest = ns0.getCpdbVersion_Dec().pyclass
getCpdbVersionResponse = ns0.getCpdbVersionResponse_Dec().pyclass
getAvailableAccessionTypesRequest = ns0.getAvailableAccessionTypes_Dec().pyclass
getAvailableAccessionTypesResponse = ns0.getAvailableAccessionTypesResponse_Dec().pyclass
mapAccessionNumbersRequest = ns0.mapAccessionNumbers_Dec().pyclass
mapAccessionNumbersResponse = ns0.mapAccessionNumbersResponse_Dec().pyclass
getAvailableFsetTypesRequest = ns0.getAvailableFsetTypes_Dec().pyclass
getAvailableFsetTypesResponse = ns0.getAvailableFsetTypesResponse_Dec().pyclass
getDefaultBackgroundSizeRequest = ns0.getDefaultBackgroundSize_Dec().pyclass
getDefaultBackgroundSizeResponse = ns0.getDefaultBackgroundSizeResponse_Dec().pyclass
overRepresentationAnalysisRequest = ns0.overRepresentationAnalysis_Dec().pyclass
overRepresentationAnalysisResponse = ns0.overRepresentationAnalysisResponse_Dec().pyclass
enrichmentAnalysisRequest = ns0.enrichmentAnalysis_Dec().pyclass
enrichmentAnalysisResponse = ns0.enrichmentAnalysisResponse_Dec().pyclass
getCpdbIdsInFsetRequest = ns0.getCpdbIdsInFset_Dec().pyclass
getCpdbIdsInFsetResponse = ns0.getCpdbIdsInFsetResponse_Dec().pyclass


if __name__ == "__main__":
    main()
