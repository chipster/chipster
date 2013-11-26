# TOOL input-output.py: InputOutputPython (Empty analysis for testing Python support..)
# INPUT input.tsv: input.tsv TYPE GENERIC 

print 'Works!'


# OUTPUT output.tsv: output.tsv 
import csv

with open('input.tsv', 'rb') as infile:
	reader = csv.reader(infile, delimiter='\t')
	with open('output.tsv', 'wb') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		for row in reader:
			writer.writerow(row)