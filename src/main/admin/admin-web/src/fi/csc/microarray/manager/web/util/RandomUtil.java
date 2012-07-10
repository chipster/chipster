package fi.csc.microarray.manager.web.util;

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Random;

public class RandomUtil {



	public static Date getRandomDate(Random rnd, int startYear) {
		Calendar cal = new GregorianCalendar();

		cal.set(rnd.nextInt(2012 - startYear) + startYear, rnd.nextInt(12) + 1, rnd.nextInt(28) + 1, rnd.nextInt(24), rnd.nextInt(60), rnd.nextInt(60));

		return cal.getTime();
	}
	
	public static Date getRandomDateToday(Random rnd) {
		Calendar cal = new GregorianCalendar();

		cal.set(cal.get(Calendar.YEAR), cal.get(Calendar.MONTH), cal.get(Calendar.DAY_OF_MONTH), rnd.nextInt(24), rnd.nextInt(60));

		return cal.getTime();
	}


	private static final String[] usernames = new String[] { "korhone", "virtane", "makine", "niemine", "makela", "hamalai", "laine", "heikkin", "koskine", "jarvine" };
	private static final String[] character = new String[] { "a", "e", "h", "i", "j", "k", "l", "m", "n", "o" };

	public static String getRandomUserName(Random rnd) {

		return usernames[rnd.nextInt(usernames.length)] + character[rnd.nextInt(character.length)];
	}

	private static final String[] sessionParts = new String[] { "hg19", "S2", "CSA", "chr", "_", "-", "gene", "genomes", "seq", "SEQ", "genomes", "analysis", "worflow", "all", "analysis", "align", "rat", "_vs_", "3", "6", "8", "2", "7" };
	
	public static String getRandomSessionName(Random rnd) {

		int count = rnd.nextInt(3) + 3;
		
		String name = "";
		for (int i = 0; i < count; i++) {
			name += sessionParts[rnd.nextInt(sessionParts.length)];
		}
		return name;
	}
	
	private static final String[] compHosts = new String[] { "hippu1.csc.fi", "hippu2.csc.fi" };
	
	public static String getRandomComp(Random rnd) {

		return compHosts[rnd.nextInt(compHosts.length)];
	}


	private static final String[] operations = new String[] { 
		"acgh-plot-combined-expression.R", 
		"classification-knn.R", 
		"norm-cdna.R", 
		"qc-agilent.R", 
		"acgh-plot-profile.R", 
		"acgh-count-overlapping-cnvs.R", 
		"stat-two-groups.R", 
		"ngs-find-nearest-genes.R", 
		"bwasw.R", 
		"FindOverlappingTool.java", 
		"stat-ROTS.R", 
		"qc-illumina.R", 
		"acgh-expression-test.R", 
		"stat-hyperG-cytoband.R", 
		"stat-hyperG-GO.R", 
		"cna-pool-bins.R", 
		"export-tab2mage.R", 
		"norm-illumina.R", 
		"bowtie-paired-end.R", 
		"ngs-annotate-miRNA-targets.R", 
		"norm-illumina-SNP.R", 
		"qc-affy-rle-nuse.R", 
		"dea-cufflinks.R", 
		"fastx-fastq-to-fasta.R", 
		"stat-hyperG-KEGG-PFAM.R", 
		"convert-bam-to-edger.R", 
		"plot-boxplot.R", 
		"generate-phenodata.R", 
		"acgh-identify-regions.R", 
		"calculate-descriptive-statistics.R", 
		"stat-linear-modelling.R", 
		"samtools-sort-index-BAM.R", 
		"norm-gene-average.R", 
		"acgh-add-cytobands.R", 
		"sort-samples.R", 
		"qc-affy-exon-rle-nuse.R", 
		"stat-one-group.R", 
		"import-ArrayExpress.R", 
		"classification.R", 
		"estimate-samplesize.R", 
		"up-down-analysis-mirna.R", 
		"fastx-collapser.R", 
		"bedtools-overlap.R", 
		"norm-affy-exon.R", 
		"MEDIPS_two_samples.R", 
		"bedtools-closestbed.R", 
		"norm-prenormalized-affy.R", 
		"plot-heatmap.R", 
		"filter-by-column-term.R", 
		"extract-features-from-BED.R", 
		"convert-RNA-bam-to-edger.R", 
		"change-interpretation.R", 
		"annotate-genelist2html.R", 
		"bedtools-bamtobed.R", 
		"ngs-dea-edger-RNA.R", 
		"ngs-filter-results-column.R", 
		"promoter-retrprom.R", 
		"samtools-idxstats.R", 
		"fastx-quality-filter.R", 
		"cna-define-experiment.R", 
		"htseq-count-own-gtf.R", 
		"acgh-call-aberrations.R", 
		"norm-specific-samples.R", 
		"samtools-snp-indel-multiple.R", 
		"ngs-filter-annotations.R", 
		"ngs-find-unique-genes.R", 
		"convert-miRBase-bam-to-edger.R", 
		"prinseq-N-filter.R", 
		"stat-adjust-pvalues.R", 
		"search-queryword.R", 
		"test-data-in.R", 
		"stat-chisq-snp.R", 
		"prinseq-AT-trimmer.R", 
		"norm-lme.R", 
		"bedtools-intersectbed.R", 
		"acgh-cluster.R", 
		"ngs-find-peaks-macs-one.R", 
		"bedtools-mergebed.R", 
		"pathways-mirna-hyperg-kegg.R", 
		"ngs-pathways-mirna-hyperg-go.R", 
		"bedtools-coveragebed.R", 
		"samtools-merge.R", 
		"mafft-einsi.sadl", 
		"filter-expression.R", 
		"ngs-find-motifs-jaspar.R", 
		"extract-genes-from-go.R", 
		"norm-chip-average.R", 
		"CombineRegionsTool.java", 
	"ngs-dea-edger.R" };

	public static String getRandomOperation(Random rnd) {

		return operations[rnd.nextInt(operations.length)];
	}

	private static final String output = "> column <- \"group\"\n" + 
			"> normalization <- \"yes\"\n" + 
			"> replicates <- \"no\"\n" + 
			"> fitting_method <- \"fit-only\"\n" + 
			"> dispersion_estimate <- \"local\"\n" + 
			"> p.value.adjustment.method <- \"BH\"\n" + 
			"> p.value.cutoff <- 0.05\n" + 
			"> image_width <- 600\n" + 
			"> image_height <- 600\n" + 
			"> # TOOL dea-deseq-RNA.R: \"Differential expression analysis using DESeq\" (This tool will perform an analysis for differentially expressed sequences using the R implementation of the DESeq algorithm.)\n" + 
			"> # INPUT data.tsv TYPE GENERIC\n" + 
			"> # INPUT phenodata.tsv TYPE GENERIC\n" + 
			"> # OUTPUT OPTIONAL de-list.tsv\n" + 
			"> # OUTPUT OPTIONAL de-list.bed\n" + 
			"> # OUTPUT OPTIONAL ma-plot.pdf\n" + 
			"> # OUTPUT OPTIONAL dispersion-plot.pdf\n" + 
			"> # OUTPUT OPTIONAL p-value-plot.pdf\n" + 
			"> # PARAMETER column: \"Column describing groups\" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)\n" + 
			"> # PARAMETER normalization: \"Apply normalization\" TYPE [yes, no] DEFAULT yes (If enabled, a normalization factor based on estimated library size is calculated.)\n" + 
			"> # PARAMETER replicates: \"Disregard replicates\" TYPE [yes, no] DEFAULT no (In order to estimate the biological and experimental variability of the data in one experiment it is necessary to have independent biological replicates of each experiment condition. However, for various reasons, biological replicates may be available for only one of the conditions or not available at all. In the former scenario, DESeq will estimate variability using the replicates of the single condition for which they are available. It is important to note that this is only an approximation and the reliability of results may suffer as a consequence. In the case where there are no replicates at all the variance is estimated by assuming the single samples from the different conditions to be replicates. The approximation will be even less reliable and results affected accordingly.)\n" + 
			"> # PARAMETER fitting_method: \"Dispersion method\" TYPE [maximum: \"fit all\", fit-only: \"fit low\"] DEFAULT maximum (The dispersion of counts for any given sequence can either be replaced with the fitted value from the dispersion model or replaced only if the fitted value is larger than the original dispersion estimate, which is the default option. The latter option optimises the balance between false positives and false negatives whereas the former minimises false positives and is therefore more conservative.)\n" + 
			"> # PARAMETER dispersion_estimate:\"Dispersion estimate\" TYPE [parametric: \"parametric\", local: \"local\"] DEFAULT local (The dispersion can be estimated either using a local fit, which is suitable in most cases - including when there are no biological independent replicate samples - or using a two-coefficient parametric model, which may be preferable under certain circumstances.)\n" + 
			"> # PARAMETER p.value.adjustment.method: \"Multiple testing correction\" TYPE [none, bonferroni: \"Bonferroni\", holm: \"Holm\", hochberg: \"Hochberg\", BH: \"BH\", BY: \"BY\"] DEFAULT BH (Multiple testing correction method.)\n" + 
			"> # PARAMETER p.value.cutoff: \"P-value cutoff\" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for statistical significance.)\n" + 
			"> # PARAMETER image_width: \"Plot width\" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)\n" + 
			"> # PARAMETER image_height: \"Plot height\" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)\n" + 
			">\n" + 
			">\n" + 
			"> ############################################################\n" + 
			"> #                                                          #\n" + 
			"> # Analaysis workflow using DESeq for normalization and     #\n" + 
			"> # statistical testing for finding differentially expressed #\n" + 
			"> # sequence tags                                            #\n" + 
			"> #                                                          #\n" + 
			"> # MG, 7.2.2012                                             #\n" + 
			"> #                                                          #\n" + 
			"> ############################################################\n" + 
			">\n" + 
			"> # Loads the libraries\n" + 
			"> library(DESeq)\n" + 
			"Loading required package: Biobase\n" + 
			"\n" + 
			"Welcome to Bioconductor\n" + 
			"\n" + 
			"  Vignettes contain introductory material. To view, type\n" + 
			"  \'browseVignettes()\'. To cite Bioconductor, see\n" + 
			"  \'citation(\"Biobase\")\' and for packages \'citation(\"pkgname\")\'.\n" + 
			"\n" + 
			"Loading required package: locfit\n" + 
			"Loading required package: akima\n" + 
			"Loading required package: lattice\n" + 
			"locfit 1.5-6 2010-01-20\n" + 
			">\n" + 
			"> # Set parameters for testing\n" + 
			"> # column <- \"group\"\n" + 
			"> # replicates <- \"yes\"\n" + 
			"> # normalization <- \"yes\"\n" + 
			"> # fitting_method <- \"maximum\"\n" + 
			"> # dispersion_estimate <- \"parametric\"\n" + 
			"> # p.value.adjustment.method <- \"BH\"\n" + 
			"> # p.value.cutoff <- 0.1\n" + 
			"> # image_height <- 600\n" + 
			"> # image_width <- 600\n" + 
			">\n" + 
			"> # Loads the normalized data\n" + 
			"> file <- c(\"data.tsv\")\n" + 
			"> dat <- read.table(file, header=T, sep=\"\t\", row.names=1)\n" + 
			">\n" + 
			"> # Separates expression values and flags\n" + 
			"> annotations <- dat[,-grep(\"chip\", names(dat))]\n" + 
			"> dat2 <- dat[,grep(\"chip\", names(dat))]\n" + 
			">\n" + 
			"> # Test needs a parameter \"groups\" that specifies the grouping of the samples\n" + 
			"> phenodata <- read.table(\"phenodata.tsv\", header=T, sep=\"\t\")\n" + 
			"> groups <- as.character (phenodata[,pmatch(column,colnames(phenodata))])\n" + 
			"> group_levels <- levels(as.factor(groups))\n" + 
			"> number_samples <- length(groups)\n" + 
			">\n" + 
			"> # If the library_size column contains data then use that as estimate\n" + 
			"> lib_size <- as.numeric(phenodata$library_size)\n" + 
			"> if (is.na(lib_size[1])) estimate_lib_size <- \"TRUE\" else estimate_lib_size <- \"FALSE\"\n" + 
			"> lib_size <- lib_size/mean(lib_size)\n" + 
			">\n" + 
			">\n" + 
			"> # Sanity checks\n" + 
			"> # only 2 group comparison is supported\n" + 
			"> if (length(unique(groups))==1 | length(unique(groups))>=3) {\n" + 
			"+ stop(\"CHIPSTER-NOTE: You need to have exactly two groups to run this analysis\")\n" + 
			"+ }\n" + 
			"> # if no biological replicates, force blind mode in dispersion estimation\n" + 
			"> if (number_samples == 2 && replicates == \"no\")  {\n" + 
			"+ stop(\"CHIPSTER-NOTE: You need to have independent biological replicates for at least one of the experiment conditions in order to reliably estimate dispersion. Alternatively, run the analysis with the disregard replicates parameter set to yes, but note that statistical power may be significantly reduced and the false positive rate may increase.\")\n" + 
			"+ }\n" + 
			"> if (number_samples == 2 && replicates == \"yes\")  {\n" + 
			"+ blind_dispersion <- TRUE\n" + 
			"+ } else {\n" + 
			"+ blind_dispersion <- FALSE\n" + 
			"+ }\n" + 
			">\n" + 
			"> # Create a counts data object\n" + 
			"> counts_data <- newCountDataSet(dat2, groups)\n" + 
			">\n" + 
			"> # Calculate scaling factors based on estimated library size, unless it is given in phenodata\n" + 
			"> # If normalization is turned off, set the size factors to 1 for all samples\n" + 
			"> if (normalization == \"yes\") {\n" + 
			"+ if (estimate_lib_size) {\n" + 
			"+ counts_data <- estimateSizeFactors(counts_data)\n" + 
			"+ } else {\n" + 
			"+ counts_data <- estimateSizeFactors(counts_data)\n" + 
			"+ sizeFactors(counts_data) <- lib_size\n" + 
			"+ }\n" + 
			"+ } else {\n" + 
			"+ sizeFactors(counts_data) <- 1\n" + 
			"+ }\n" + 
			">\n" + 
			"> # Estimate dispersion values for each gene and replace with fitted values\n" + 
			"> # Use sharingMode parameter to control how conservative the replacement will be\n" + 
			"> # Use fitType to control for parametric or local fit\n" + 
			"> if (blind_dispersion) {\n" + 
			"+ counts_data <- estimateDispersions(counts_data, method=\"blind\", sharingMode=\"fit-only\",\n" + 
			"+ fitType=dispersion_estimate)\n" + 
			"+ } else {\n" + 
			"+ counts_data <- estimateDispersions(counts_data, method=\"pooled\", sharingMode=fitting_method,\n" + 
			"+ fitType=dispersion_estimate)\n" + 
			"+ }\n" + 
			">\n" + 
			">\n" + 
			"> # Function that produces a qc plot to check dispersion estimates\n" + 
			"> plotDispEsts <- function(cds) {\n" + 
			"+ plot(rowMeans( counts(cds, normalized=TRUE)), fitInfo(cds)$perGeneDispEsts,pch = \'.\',\n" + 
			"+ log=\"xy\", main=\"Dispersion plot\", xlab=\"normalized counts\", ylab=\"dispersion\")\n" + 
			"+ xg <- 10^seq( -.5, 5, length.out=300)\n" + 
			"+ lines(xg, fitInfo(cds)$dispFun(xg), col=\"red\")\n" + 
			"+ legend(x=\"topright\", legend=\"fitted dipersion\", col=\"red\", cex=1, pch=\"-\")\n" + 
			"+ }\n" + 
			">\n" + 
			"> # Make dispersion plot\n" + 
			"> pdf(file=\"dispersion-plot.pdf\")\n" + 
			"> plotDispEsts(counts_data)\n" + 
			"Warning message:\n" + 
			"In xy.coords(x, y, xlabel, ylabel, log) :\n" + 
			"  860 y values <= 0 omitted from logarithmic plot\n" + 
			"> dev.off()\n" + 
			"null device\n" + 
			"          1\n" + 
			">\n" + 
			"> # Calculate statistic for differential expression\n" + 
			"> results_table <- nbinomTest(counts_data, group_levels[2], group_levels[1] )\n" + 
			">\n" + 
			"> # Merge with original data table\n" + 
			"> output_table <- cbind (dat2, results_table[,-1])\n" + 
			">\n" + 
			"> # Adjust p-values\n" + 
			"> output_table$padj <- p.adjust(output_table$pval, method=p.value.adjustment.method)\n" + 
			">\n" + 
			"> # Filter out the significant ones\n" + 
			"> significant_table <- output_table[ (output_table$padj <  p.value.cutoff),]\n" + 
			">\n" + 
			"> # Remove rows with NA adjusted p-values\n" + 
			"> significant_table <- significant_table[! (is.na(significant_table$padj)),]\n" + 
			">\n" + 
			"> # Order results based on raw p-values\n" + 
			"> significant_table <- significant_table[ order(significant_table$pval), ]\n" + 
			">\n" + 
			"> # Output the table\n" + 
			"> if (dim(significant_table)[1] > 0) {\n" + 
			"+ write.table(significant_table, file=\"de-list.tsv\", sep=\"\t\", row.names=T, col.names=T, quote=F)\n" + 
			"+ }\n" + 
			">\n" + 
			"> # Make histogram of p-values with overlaid significance cutoff and uniform distribution\n" + 
			"> pdf (file=\"p-value-plot.pdf\")\n" + 
			"> hist(output_table$pval, breaks=100, col=\"blue\",\n" + 
			"+ border=\"slateblue\", freq=FALSE,\n" + 
			"+ main=\"P-value distribution\", xlab=\"p-value\", ylab=\"proportion (%)\")\n" + 
			"> hist(output_table$padj, breaks=100, col=\"red\",\n" + 
			"+ border=\"slateblue\", add=TRUE, freq=FALSE)\n" + 
			"> abline(h=1, lwd=2, lty=2, col=\"black\")\n" + 
			"> abline(v=p.value.cutoff, lwd=2, lty=2, col=\"green\")\n" + 
			"> legend (x=\"topright\", legend=c(\"p-values\",\"adjusted p-values\", \"uniform distribution\", \"significance cutoff\"), col=c(\"blue\",\"red\",\"black\",\"green\"),\n" + 
			"+ cex=1, pch=15)\n" + 
			"> dev.off()\n" + 
			"null device\n" + 
			"          1\n" + 
			">\n" + 
			"> # Define function for making MA-plot of significant findings\n" + 
			"> plotDE <- function(res)\n" + 
			"+ plot(res$baseMean, res$log2FoldChange,\n" + 
			"+ log=\"x\", pch=20, cex=.25, col = ifelse( res$padj < p.value.cutoff, \"red\", \"black\"),\n" + 
			"+ main=\"MA plot\", xlab=\"mean counts\", ylab=\"log2(fold change)\")\n" + 
			">\n" + 
			"> # Make MA-plot\n" + 
			"> pdf(file=\"ma-plot.pdf\")\n" + 
			"> plotDE(results_table)\n" + 
			"> legend (x=\"topleft\", legend=c(\"significant features\",\"not significant\"), col=c(\"red\",\"black\"),\n" + 
			"+ cex=1, pch=19)\n" + 
			"> dev.off()\n" + 
			"null device\n" + 
			"          1\n" + 
			">\n" + 
			"> # EOF";


	public static String getRandomOutputtext(Random rnd) {
		int start = rnd.nextInt(output.length());
		int end = rnd.nextInt(output.length());

		return output.substring(Math.min(start, end), Math.max(start, end));
	}

	private static final String[] errors = new String[] { 

		"in getBM(filters = \"ensembl_gene_id\", values = ensembl_id_list,  :\n" + 
				"Invalid attribute(s): hgnc_symbol\n" + 
				"Please use the function \'listAttributes\' to get valid attribute names\n" + 
				"Starting R failed.\n", 
				"inwrite.table(results_list, file = \"de-genes-cufflinks.tsv\", sep = \"\t\",  :\n" + 
						"invalid \'row.names\' specification\n", 
						"inread.table(file, header = F, sep = \"\t\") :\n" + 
								"no lines available in input\n", 
								"insort.list(y) :\n" + 
										"\'x\' must be atomic for \'sort.list\'\n" + 
										"Have you called \'sort\' on a list?\n" + 
										"Calls: factor -> sort.list\n", 
										"inboxplot.default(as.data.frame(dat2), las = 2, col = sample_colots) :\n" + 
												"object \'sample_colots\' not found\n" + 
												"Calls: boxplot -> boxplot.default\n", 
												"inif (file.type == \"application/x-gzip\") { :\n" + 
														"argument is of length zero\n" + 
														"Calls: unzipIfGZipFile -> isGZipFile\n", 
														"inaddition: Warning message:\n" + 
																"running command \'file -Lb --mime-type fastqfile\' had status 1\n", 
																"incmdscale(as.dist(dd), k = ndim) :\n" + 
																		"\'k\' must be in {1, 2, ..  n - 1}\n" + 
																		"Calls: plotMDS.dge -> cmdscale\n" + 
																		"Running R script failed.\n", 
																		"ingetGoToEntrezMap_db(p) :\n" + 
																				"The genes you are testing do not have any corresponding GO terms for the ontology you are searching.\n" + 
																				"Calls: hyperGTest ... categoryToEntrezBuilder -> categoryToEntrezBuilder -> getGoToEntrezMap_db\n" + 
																				"Building the peak model failed. Retry by lowering the m-fold value or enabling the automatic m-fold adjustment.\n", 
																				"incat(\"Out of the\", number_genes_tested, \"genes tested there were no statistically significantly differentially expressed ones found.\") :\n" + 
																						"object \'number_genes_tested\' not found\n" + 
																						"It appears as though the input file(s) are not in fastq format. Please check input files or rerun the tool but with the \'Input file format\' parameter set to \'FASTA\'.\n", 
																						"innames(x) <- value :\n" + 
																								"\'names\' attribute [13] must be the same length as the vector [4]\n" + 
																								"Calls: colnames<-\n" + 
																								": You need to have exactly two groups to run this analysis\n", 
																								"inqtclust(dat3, radius = rad) :\n" + 
																										"All points in one cluster, try smaller radius.\n", 
																										"inlm.fit(design, t(M)) : incompatible dimensions\n" + 
																												"Calls: lmFit -> lm.series -> lm.fit\n" + 
																												"It appears that the two input files are not matepairs. Please checkthat the correct input files were selected.\n", 
																												"inread.table(file = \"bam-out.tsv\", header = FALSE, sep = \"\t\") :\n" + 
																														"no lines available in input\n" + 
																														"Required output file is missing.\n", 
																														"infile(file, \"r\", encoding = encoding) :\n" + 
																																"cannot open the connection\n" + 
																																"Calls: source -> file\n", 
																																"inaddition: Warning message:\n", 
																																"infile(file, \"r\", encoding = encoding) :\n" + 
																																		"cannot open file \'/opt/chipster/comp/modules/common/R-2.12/zip-utils.R\': No such file or directory\n", 
																																		"inqtclust(dat3, radius = rad) :\n" + 
																																				"Could not find a valid clustering, try again with different radius.\n", 
																																				"inparametricDispersionFit(means, disps) :\n" + 
																																						"Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See \'?estimateDispersions\')\n" + 
																																						"Calls: estimateDispersions ... estimateAndFitDispersionsFromBaseMeansAndVariances -> parametricDispersionFit\n", 
																																						"in`row.names<-.data.frame`(`*tmp*`, value = c(1L, 0L)) :\n" + 
																																								"invalid \'row.names\' length\n" + 
																																								"Calls: rownames<- -> row.names<- -> row.names<-.data.frame\n", 
																																								"in.checkKeys(value, Lkeys(x), x@ifnotfound) :\n" + 
																																										"value for \"lysA\" not found\n" + 
	"Calls: unlist ... keys<- -> keys<- -> Lkeys<- -> Lkeys<- -> .checkKeys" };


	public static String getRandomError(Random rnd) {

		return errors[rnd.nextInt(errors.length)];
	}
}
