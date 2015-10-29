# TOOL vep.py: "Ensembl Variant Effect Predictor" (Given a VCF file, this tool determines the effect of SNPs, insertions, deletions, CNVs or structural variants on genes, transcripts, and protein sequence, as well as regulatory regions. It uses the Ensembl VEP service running at the EBI. It works currently with human hg38 data, but more genomes will be added.)
# INPUT input_file: "VCF file" TYPE GENERIC (Input variants in VCF format)
# OUTPUT vep.tsv: "VEP file"

import urllib
import urllib2
import json
import time
import os
#from tool_utils import *

"""
Should produce the same results with the web interface, when 
"Transcript support level" and "1000 Genomes global minor allele frequency" are disabled, because the REST API doesn't give this information.

The VEP output format is really complicated. It might be a good idea to use something else.

Order may be different, so use this command to check if results are the same:

diff <(cat web.vep | sort ) <(cat output_file | sort )
"""

def get(dict, key):
    if key in dict:
        return str(dict.get(key))
    return '-'

def append(list, prefix, dict, key):
    if key in dict:
        return list.append(prefix + str(dict.get(key)))

def writeToFile(input_region, consequence_type, feature_type, f):
    line_count = 0
    consequences = input_region.get(consequence_type)
    if not consequences:
        return 0
    for conseq in consequences:
        vep_line = []
        uploaded_variation = '.'
        location = input_region['seq_region_name'] + ':' + str(input_region['start']) + '-' + str(input_region['end'])
        allele = get(conseq, 'variant_allele')
        consequence = ','.join(conseq['consequence_terms'])
        impact = conseq['impact']
        symbol = get(conseq, 'gene_symbol')
        gene = get(conseq, 'gene_id')
        feature = get(conseq, 'transcript_id')
        if 'regulatory_feature_id' in conseq:
            feature = get(conseq, 'regulatory_feature_id')
        biotype = get(conseq, 'biotype')
        exon = get(conseq, 'exon')
        intron = get(conseq, 'intron')
        hgvsc = '-'
        hgvsp = '-'
        cdna_position = get(conseq, 'cdna_start')
        cds_position = get(conseq, 'cds_start')
        protein_position = get(conseq, 'protein_start')
        amino_acids = get(conseq, 'amino_acids')
        codons = get(conseq, 'codons')
    
        extras = []            
        append(extras, 'DISTANCE=', conseq, 'distance')
        append(extras, 'STRAND=', conseq, 'strand')
        append(extras, 'SYMBOL_SOURCE=', conseq, 'gene_symbol_source')
        append(extras, 'HGNC_ID=', conseq, 'hgnc_id')
        if 'sift_prediction' in conseq:
            extras.append('SIFT=' + conseq['sift_prediction'] + '(' + str(conseq['sift_score']) + ')')
        if 'polyphen_prediction' in conseq:
            extras.append('PolyPhen=' + conseq['polyphen_prediction'] + '(' + str(conseq['polyphen_score']) + ')')

        existing_variation = '-'
        variants = []
        somatics = []
        phenos = []
        colocated = input_region.get('colocated_variants')
        if colocated:
            for variant in colocated:
                variants.append(variant['id'])                  
                append(somatics, '', variant, 'somatic')
                append(phenos, '', variant, 'phenotype_or_disease')
            if variants:
                existing_variation = ','.join(variants)
            if somatics:
                extras.append('SOMATIC=' + ','.join(somatics))
            if phenos:
                extras.append('PHENO=' + 'j'.join(phenos))

        extra = ';'.join(extras)

        vep_line = [uploaded_variation, location, allele, consequence, impact, symbol, gene, feature_type, feature, biotype, exon, intron, hgvsc, hgvsp, cdna_position, cds_position, protein_position, amino_acids, codons, existing_variation, extra]

        f.write('\t'.join(vep_line) + '\r\n')
        line_count += 1
    return line_count

# parse VCF file and create an array of variants
def parse_lines(lines):
    variants = []
    for line in lines:
        if line.startswith('#'):
            continue
        columns = line.split()
        chrom = columns[0]
        pos = columns[1]
        ref = columns[3]
        alt = columns[4]
        end_pos = str(int(pos) + len(ref) - 1)
        variants.append(chrom + ' ' + pos + ' ' + end_pos + ' ' + ref + ' ' + alt)

    return variants

def query(lines):
    url = 'http://rest.ensembl.org/vep/homo_sapiens/region?numbers=true'

    variants = parse_lines(lines)
    # input format example: { "variants" : ["21 26960070 rs116645811 G A . . .", "21 26965148 rs1135638 G A . . ." ] }
    values = {'variants': variants }
    # convert to json string  
    data = json.dumps(values)
    # send request
    req = urllib2.Request(url, data)
    req.add_header("Content-type", "application/json")
    response = urllib2.urlopen(req)    
    # parse response json
    return json.loads(response.read())

def write(response_data):
    line_count = 0
    with open('vep.tsv', 'a') as f:
        # write file header
        f.write('#Uploaded_variation	Location	Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra\r\n')

        # write consequences of each input region
        for input_region in response_data:
            # there are many types of consequences
            line_count += writeToFile(input_region, 'transcript_consequences', 'Transcript', f)
            line_count += writeToFile(input_region, 'intergenic_consequences', '-', f)
            line_count += writeToFile(input_region, 'motif_feature_consequences', '-', f)
            line_count += writeToFile(input_region, 'regulatory_feature_consequences', 'RegulatoryFeature', f)
    print('wrote ' + str(line_count) + ' output lines')

def main():

    # the file may exist from previous runs, when debugging this locally
    if os.path.exists('vep.tsv'):
        os.remove('vep.tsv')
    
    lines = []
    with open('input_file') as f:
        for line in f:
            lines.append(line)
            if len(lines) >= 1000:
                print('query ' + str(len(lines)) + ' variants')
                write(query(lines))
                lines = []
                # REST API allows 15 requests per second, so this sleep makes sure 
                # that we can safely run about 30 instances of this tool in parallel without having
                # to wait for the quota reset (which may take up to an hour).
                # https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
                time.sleep(2)      
        print('query ' + str(len(lines)) + ' variants')
        write(query(lines))      

"""
Transform input JSON: 

[
  {
    "input": "21  26960070  rs116645811 G A . . .",
    "assembly_name": "GRCh38",
    "end": 26960070,
    "seq_region_name": "21",
    "transcript_consequences": [
      {
        "gene_id": "ENSG00000154736",
        "variant_allele": "A",
        "biotype": "protein_coding",
        "gene_symbol_source": "HGNC",
        "consequence_terms": [
          "intron_variant"
        ],
        "strand": -1,
        "hgnc_id": "HGNC:221",
        "gene_symbol": "ADAMTS5",
        "transcript_id": "ENST00000284987",
        "impact": "MODIFIER"
      }
    ],
    "strand": 1,
    "id": "rs116645811",
    "allele_string": "G/A",
    "most_severe_consequence": "intron_variant",
    "start": 26960070
  }
]

To VEP format:

.	1:142492-142492	A	upstream_gene_variant	MODIFIER	RP11-34P13.16	ENSG00000269981	Transcript	ENST00000595919	processed_pseudogene	-	-	-	-	-	-	-	-	-	rs202086374	DISTANCE=4527;STRAND=-1;SYMBOL_SOURCE=Clone_based_vega_gene
"""

if __name__ == "__main__":
    main()

