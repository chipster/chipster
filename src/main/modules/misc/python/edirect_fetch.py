# TOOL edirect_fetch.py: "Retrieve sequences from NCBI" (Retrieves sequences from NCBI databases using different search criteria.)
# OUTPUT OPTIONAL sequences.txt: sequences.txt (Retrieved sequence set)
# OUTPUT OPTIONAL sequences.fasta: sequeces.fasta (Result sequences in fasta format)
# OUTPUT OPTIONAL edirect.log
# PARAMETER db: "Sequence type" TYPE [protein: Protein, nucleotide: Nucleotide] DEFAULT protein
# PARAMETER q1field: "Search field for first query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER q1term: "Query term" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL log_op: "Logical operator" TYPE [AND: AND, OR: OR, NOT: NOT] DEFAULT AND (Logical operators used to combine the search terms)
# PARAMETER OPTIONAL q2field: "Search field for second query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name", SLEN: "Sequence length"] DEFAULT ALL (Select the search field )
# PARAMETER OPTIONAL q2term: "Second search term" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL q3field: "Search field for third query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name", SLEN: "Sequence length"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q3term: "Third search term" TYPE STRING (Search term or word.)
# PARAMETER outformat: "Output format" TYPE [fasta: FASTA, gp: "Genbank proteins", gb: "Genbank"] DEFAULT fasta (Logical operators used to build the search terms)


#To make edirect work:
# sudo aptitude install libhttp-parser-perl
#

import socket
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET

def main():
       
    edirect_path = "/opt/chipster/tools/edirect/"
    esearch_path = edirect_path + "esearch"
    efetch_path = edirect_path + "efetch"

    if outformat == "fasta":
        outfile = "sequences.fasta"
    else:
        outfile = "sequences.txt"
 
    # each argument must be a separate list item
    query = [esearch_path, "-db", db, "-query"]

    # construct a query string
    query_str = q1term + "[" + q1field + "]"
    if q2term:
        query_str += " " + log_op + " " + q2term + "[" + q2field + "]"
    if q3term:
        query_str += " " + log_op + " " + q3term + "[" + q3field + "]"

    # append query string to the argument list
    query.append(query_str)
    
    # run the query command
    # save standard output and error streams to variables xml and err
    esearch_process = Popen(query, stdout=PIPE, stderr=PIPE)
    xml, err = esearch_process.communicate()
    if err:
        # something in error stream, give up
        raise Exception("Error in esearch: " + err)

    # parse xml
    str_hits = ET.fromstring(xml).find("Count").text    
    num_hits = int(str_hits)

    if num_hits == 0:
        raise Exception("CHIPSTER-NOTE: No hits found for query: esearch " + " ".join(query))

    elif num_hits > 50000:
        raise Exception("CHIPSTER-NOTE: Query produced more than 50000 hits. Number of hits " + str_hits)
        
    # fetch sequences
    # read standard input from a variable xml and write output to a file, because it may be big
    with open("hits.txt", "w") as hits:
        efetch_process = Popen([efetch_path, "-format", outformat], stdout=hits, stdin=PIPE)
        efetch_process.communicate(input=xml)

    # read the result file, filter it and write to another file
    with open("hits.txt", "r") as hits, open(outfile, "w") as out:
        for line in hits:
            # remove empty lines
            if line != "\n":                
                out.write(line)

if __name__ == "__main__":
    main()
