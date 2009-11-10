#################################################
# 
# Readme for R scripts in Chipster version 1.3.0
# Jarno Tuimala, 29th January 2009
#
#################################################


#################
# 
# New tools
#
#################

- Add basic annotations to data (annotate-add-to-data.R) adds the annotations to the end of the
  data file. This enables one to filter the data file using the Filter by column -tool. This way
  it is possible to select genes, e.g., only from a certain KEGG pathway, GO ontology class or
  chromosomal cytoband.
- Hypergeometric test for cytobands (stat-hyperG-cytoband.R) runs the hypergeometric test for
  gene enrichment for certain chromosomal areas. 
- Classification (classification.R) enables one to run a class prediction / supervised clustering /
  classification analysis using all sorts of methods. This tool makes the older tool KNN classification
  replete, except that the newest tool does not implement the validation on test set yet.
- Process pre-normalized (norm-prenormalized.R) converts the prenormalized data set that has been
  imported through import tool (if it is not in Chipster-normalized format) into Chipster-normalized
  format, most importantly, adds text "chip." in front of the expression value columns, and generates
  an empty phenodata object.
- Sort genes (sort-genes.R) sort the data rowwise in ascending to descending order according to some
  selected data column. This can be used in conjunction with extract rows -tool to select a specified
  number of genes from the beginning or end of the data.
- Extract genes (extract-genes.R) extract a specified number of genes from the top of the data table.
  Works especially well after the data has been sorted using some column. One example of a workflow
  could be to first calculate descriptive statistics, sort the data using the standard deviation column,
  and then extract the 50 top-most genes.
- Normalize to specific genes (norm-specific-genes.R) allows the user to normalize to specific (positive)
  control genes. An average of the specified genes on each chip is calculated and subtracted from all 
  the values on the same chip. It requires two inputs: the normalized data and a list of gene identifiers
  where the first column is titles identifier, and hold the row names used in the data. The second column
  should be titled chip, and should be empty.


#################
# 
# Recent changes
#
#################

- Hierarchical clustering (under clustering category) was using uncentered Pearson correlation
  instead of the intended centered Pearson correlation. This has been corrected, and the tool
  should now as originally intended. This also makes the bootstrapping result, and the tree
  more comparable, because they are generated using the same distance method.
- Filter by column -tool can now be used for selecting rows (genes) using some key words, also.
  Thanks for Amanda Miotto, Australia, for the modification!
- Normalize cDNA -tools sometimes resulted in replacing gene identifiers with a running number.
  This is now hopefully fixed. The same phenomenan still happens, if the Agilent tools are applied
  to generic cDNA data, but it might not be worth correcting at this point.
- Fast t-test has been added to the two-sample tests. This test is equivalent to a linear regression
  analysis using one factor with two levels. The difference to the option t-test is that fast-t-test
  is calculated on the transposed data matrix using the commnd lm() from R, which makes the calculation
  lightning fast (about a thousand times faster then the standard t-test). Execution time is about the 
  same for the empiricalBayes and fast-t-test. The reason for leaving the standard t-test in Chipster is
  that the p-values given by the fast-t-test are slightly different from the ones reported by the standard
  t-test for a small number of genes.  
- Tool Calculate descriptive statistics does not expect to get a phenodata anymore. In other words,
  the tool now works for plain data tables, even without phenodata.
- Tool Filter by column does not expect to get a phenodata anymore. In other words, the tool now works 
  for plain data tables, even without phenodata.
- Adjust P-values tool now offers Storey's Q-value correction as a new option. 
- Visualization tools Correlogram, Boxplot, Heatmap, Histogram, Dendrogram now use description column from 
  phenodata for labeling samples in the plots. This also means that these tool now require phenodata as an
  input.
- K-means - test K -tool now uses 10 random starts (instead of 1 random start in the previous version).
  This might not still be enough for getting an optimal solution for a large number of clusters, but should
  give a more reasonable estimate, and run fast.
- K-means clustering -tool now uses either the number of genes or ten, depending which is smaller, as
  the number of random starts in the analysis. Previously just one random start was used.
- Calculate descriptive statistics crashed on execution on chips. If chips was chosen, no file descriptives.tsv
  was generated, and that made the tool to fail. Should be corrected now.
- Normalize Affymetrix exon -tool contained wrong CDF environment names, when gene summaries were requested. 
  Corrected to work with Entrez Gene annotations.
- Change interpretations now supports tranformation to and forth to log10 and ln and linear data. The tools
  does not require phenodata to run anymore.


#################
# 
# Known bugs
#
#################

- Flags are not returned for Agilent, cDNA or Illumina chips, even if specified. 
- Combine probes to genes -tool crashes when generating identifiers for the genes.


#################
# 
# Notes for users
#
#################

- Possibility to tell how large GO categories or KEGG pathways are tested during the gene enrichment analysis
  has not been implemented, since this functionality relates to filtering the results of the analysis and
  not restricting the analysis as such. In other words, even if we added the option to restrict the minimum
  category size, it would not affect the analysis, which would be run for all categories, but only how the
  results are reported. Therefore, we thought it would better to output all significant results, and let the 
  user decide which ones are of biological interest.
- Gene enrichment analyses in Chipster do not return the genes that belong to the significant categories. 
  Implementing this would make the tool run even slower (much more so), but we have crated a new annotation tool
  that adds the annotation to the data file. That way the genes can be sorted (or filtered) in the spreadsheet view 
  according to the GO classes or KEGG pathways into which they fall.
- Instead of implementing a more complete search tool, we suggest to use the Add basic annotations to data and
  Filter by column (using a keyword) tools to perform the same operation. If you need to search for several genes
  at the same time (filter by gene list), create a simple text file: first line should read as "identifier", and
  the following ones should contain one gene identifier (the ones used as row names in Chipster) each. Import
  the file to Chipster, and put the data file and the gene list file into a Venn diagram to perform the filtering. 