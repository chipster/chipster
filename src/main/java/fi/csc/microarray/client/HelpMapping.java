package fi.csc.microarray.client;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.operation.OperationDefinition;

public class HelpMapping {

	private static final String DEFAULT_HELP_PAGE = "chipster-manual/tools.html";
	private static Map<String, String> mappings = new HashMap<String, String>();

	static {
		mappings.put("Preprocessing/Filter by CV", "chipster-manual/filter-cv.html");
		mappings.put("Preprocessing/Filter by expression", "chipster-manual/filter-expression.html");
		mappings.put("Preprocessing/Filter by flags", "chipster-manual/filter-flags.html");
		mappings.put("Preprocessing/Filter by interquartile range", "chipster-manual/filter-iqr.html");
		mappings.put("Preprocessing/Filter by standard deviation", "chipster-manual/filter-sd.html");
		mappings.put("Preprocessing/Filter using a column term", "chipster-manual/filter-by-column-term.html");
		mappings.put("Preprocessing/Filter using a column value", "chipster-manual/filter-by-column-value.html");
		mappings.put("Preprocessing/Impute missing values", "chipster-manual/impute.html");
		mappings.put("Preprocessing/Remove missing values", "chipster-manual/na-omit.html");

		mappings.put("Quality control/Affymetrix basic", "chipster-manual/qc-affy.html");
		mappings.put("Quality control/Affymetrix - using RLE and NUSE", "chipster-manual/qc-affy-rle-nuse.html");
		mappings.put("Quality control/Agilent 1-color", "chipster-manual/qc-agilent-one-color.html");
		mappings.put("Quality control/Agilent 2-color", "chipster-manual/qc-agilent.html");
		mappings.put("Quality control/cDNA", "chipster-manual/qc-cdna.html");
		mappings.put("Quality control/Illumina", "chipster-manual/qc-illumina.html");
							
		mappings.put("Normalisation/Affymetrix exon arrays", "chipster-manual/norm-affy-exon.html");
		mappings.put("Normalisation/Affymetrix", "chipster-manual/norm-affy.html");
		mappings.put("Normalisation/Affymetrix SNP arrays", "chipster-manual/norm-affy-snp.html");
		mappings.put("Normalisation/Affymetrix gene arrays", "chipster-manual/norm-affy-gene.html");
		mappings.put("Normalisation/Agilent miRNA", "chipster-manual/norm-agilent-mirna.html");		
		mappings.put("Normalisation/Illumina SNP arrays", "chipster-manual/norm-illumina-snp.html");
		mappings.put("Normalisation/Process prenormalized", "chipster-manual/norm-process-prenormalized.html");
		mappings.put("Normalisation/Process prenormalized affy", "chipster-manual/norm-prenormalized-affy.html");
		mappings.put("Normalisation/Agilent 1-color", "chipster-manual/norm-agilent-1color.html");
		mappings.put("Normalisation/Agilent 2-color", "chipster-manual/norm-agilent.html");
		mappings.put("Normalisation/cDNA", "chipster-manual/norm-cdna.html");
		mappings.put("Normalisation/Illumina", "chipster-manual/norm-illumina.html");
		mappings.put("Normalisation/Illumina - lumi pipeline", "chipster-manual/norm-illumina-lumi.html"); 
		mappings.put("Normalisation/Random effects", "chipster-manual/norm-lme.html");
		mappings.put("Normalisation/Normalize to chip average", "chipster-manual/norm-chip-average.html");
		mappings.put("Normalisation/Normalize to gene average", "chipster-manual/norm-gene-average.html");
		mappings.put("Normalisation/Normalize to specific samples", "chipster-manual/norm-specific-samples.html");
		mappings.put("Normalisation/Normalize to specific genes", "chipster-manual/norm-specific-genes.html");
		
		mappings.put("Statistics/Gene set test", "chipster-manual/stat-geneset.html");
		mappings.put("Statistics/One sample tests", "chipster-manual/stat-one-group.html");
		mappings.put("Statistics/Several groups tests", "chipster-manual/stat-several-groups.html");
		mappings.put("Statistics/Single-slide methods", "chipster-manual/stat-singleslide.html");
		mappings.put("Statistics/Time series", "chipster-manual/stat-timeseries.html");
		mappings.put("Statistics/Two groups tests", "chipster-manual/stat-two-groups.html");
		mappings.put("Statistics/ROTS", "chipster-manual/stat-ROTS.html");
		mappings.put("Statistics/NMDS", "chipster-manual/ordination-nmds.html");
		mappings.put("Statistics/PCA", "chipster-manual/ordination-pca.html");
		mappings.put("Statistics/Sample size estimation", "chipster-manual/stat-estimate-sample-size.html");
		mappings.put("Statistics/Sample size calculations with an adapted BH method", "chipster-manual/sample-size-with-bh.html");
		mappings.put("Statistics/Correlate with phenodata", "chipster-manual/stat-correlate-phenodata.html");
		mappings.put("Statistics/Correlate miRNA with target expression", "chipster-manual/correlate-mirna.html");
		mappings.put("Statistics/Linear modelling", "chipster-manual/stat-linear-modelling.html");
		mappings.put("Statistics/SAM", "chipster-manual/stat-sam.html");
		mappings.put("Statistics/Adjust p-values", "chipster-manual/stat-adjust-p-values.html");
		mappings.put("Statistics/Calculate descriptive statistics", "chipster-manual/calculate-descriptive-statistics.html");
		mappings.put("Statistics/DCA", "chipster-manual/ordination-dca.html");
		mappings.put("Statistics/Association analysis", "chipster-manual/stat-chisq-snp.html");
		mappings.put("Statistics/Up-down analysis of miRNA targets", "chipster-manual/up-down-analysis-mirna.html");
		
		mappings.put("Clustering/Hierarchical", "chipster-manual/cluster-hierarchical.html");
		mappings.put("Clustering/K-Means", "chipster-manual/cluster-kmeans.html");
		mappings.put("Clustering/KNN classification", "chipster-manual/cluster-knn-classification.html");
		mappings.put("Clustering/Quality Threshold (QT)", "chipster-manual/cluster-qt.html");
		mappings.put("Clustering/Self-organizing map (SOM)", "chipster-manual/cluster-som.html");
		mappings.put("Clustering/K-Means - estimate K", "chipster-manual/cluster-kmeans-testk.html");
		mappings.put("Clustering/Classification", "chipster-manual/cluster-classification.html");
		
		mappings.put("Pathways/Bayesian network", "chipster-manual/pathway-bayesian.html");
		mappings.put("Pathways/Boolean network", "chipster-manual/pathway-boolean-bestim.html");
		mappings.put("Pathways/Hypergeometric test for GO", "chipster-manual/pathways-hypergeometric-go.html");
		mappings.put("Pathways/Hypergeometric test for KEGG or PFAM", "chipster-manual/pathways-hypergeometric-kegg.html");
		mappings.put("Pathways/SAFE test for KEGG pathway enrichment", "chipster-manual/pathways-hypergeometric-safe.html");
		mappings.put("Pathways/Gene set test", "chipster-manual/stat-geneset.html");
		mappings.put("Pathways/Protein interactions from IntAct", "chipster-manual/pathways-intact.html");
		mappings.put("Pathways/Associations to Reactome pathways", "chipster-manual/pathways-reactome.html");
		mappings.put("Pathways/Hypergeometric test for ConsensusPathDB", "chipster-manual/pathways-hypergeometric-cpdb.html");
		mappings.put("Pathways/Hypergeometric test for cytobands", "chipster-manual/pathways-hypergeometric-cytobands.html");
		mappings.put("Pathways/KEGG enrichment for miRNA targets", "chipster-manual/pathways-hyperg-mirna-kegg.html");
		mappings.put("Pathways/GO enrichment for miRNA targets", "chipster-manual/pathways-hyperg-mirna-go.html");

		mappings.put("Visualisation/Boxplot", "chipster-manual/plot-boxplot.html");
		mappings.put("Visualisation/Chromosomal position", "chipster-manual/plot-chrom-pos.html");
		mappings.put("Visualisation/Correlogram", "chipster-manual/plot-correlogram.html");
		mappings.put("Visualisation/Dendrogram", "chipster-manual/plot-dendrogram.html");
		mappings.put("Visualisation/Heatmap", "chipster-manual/plot-heatmap.html");
		mappings.put("Visualisation/Histogram", "chipster-manual/plot-histogram.html");
		mappings.put("Visualisation/Idiogram", "chipster-manual/plot-idiogram.html");
		mappings.put("Visualisation/Venn diagram", "chipster-manual/plot-venn-diagram.html");
		mappings.put("Visualisation/Volcano plot from existing results", "chipster-manual/plot-volcano-data-exists.html");
		mappings.put("Visualisation/Volcano plot", "chipster-manual/plot-volcano.html");
		
		mappings.put("Promoter Analysis/Retrieve promoters", "chipster-manual/promoter-retrprom.html");
		mappings.put("Promoter Analysis/Weeder", "chipster-manual/promoter-tfbs.html");
		mappings.put("Promoter Analysis/ClusterBuster", "chipster-manual/promoter-cbust.html");
		mappings.put("Promoter Analysis/Cosmo", "chipster-manual/promoter-tfbs-cosmo.html");
		
		mappings.put("Annotation/Agilent, Affymetrix or Illumina genelist", "chipster-manual/annotate-genelist2html.html");
		mappings.put("Annotation/Find miRNA targets", "chipster-manual/annotate-miRNA-targets.html");
		mappings.put("Annotation/Add annotations to data", 	"chipster-manual/annotate-add-to-data.html");
		mappings.put("Annotation/Agilent miRNA", "chipster-manual/annotate-miRNA.html");
		
		mappings.put("Utilities/Word-based query", "chipster-manual/search-queryword.html");
		mappings.put("Utilities/Export GEO's SOFT format", "chipster-manual/export-soft.html");
		mappings.put("Utilities/Export tab2mage format", "chipster-manual/export-tab2mage.html");
		mappings.put("Utilities/Extract genes from clustering", "chipster-manual/extract-genes-from-clustering.html");
		mappings.put("Utilities/Extract genes using a p-value", "chipster-manual/extract-genes-from-stattest.html");
		mappings.put("Utilities/Extract samples from dataset", "chipster-manual/extract-samples-from-dataset.html");
		mappings.put("Utilities/Average replicate chips", "chipster-manual/average-replicates.html");
		mappings.put("Utilities/Calculate fold change", "chipster-manual/calculate-fold-change.html");
		mappings.put("Utilities/Generate phenodata", "chipster-manual/generate-phenodata.html");
		mappings.put("Utilities/Import from GEO", "chipster-manual/import-from-geo.html");
		mappings.put("Utilities/Merge tables", "chipster-manual/merge-tables.html");
		mappings.put("Utilities/Merge data sets", "chipster-manual/merge-datasets.html");
		mappings.put("Utilities/Search by correlation", "chipster-manual/search-correlation.html");
		mappings.put("Utilities/Search by gene name", "chipster-manual/search-queryword.html");
		mappings.put("Utilities/Merge tables", "chipster-manual/merge-tables.html");
		mappings.put("Utilities/Sort samples", "chipster-manual/sort-samples.html");
		mappings.put("Utilities/Delete columns", "chipster-manual/delete-columns.html");
		mappings.put("Utilities/Filter using a column", "chipster-manual/filter-by-column.html");
		mappings.put("Utilities/Combine probes to genes", "chipster-manual/combine-probes-to-genes.html");
		mappings.put("Utilities/Extract genes", "chipster-manual/extract-genes.html");
		mappings.put("Utilities/Change interpretation",	"chipster-manual/change-interpretation.html");
		mappings.put("Utilities/Sort genes", "chipster-manual/sort-genes.html");
		mappings.put("Utilities/Random sampling", "chipster-manual/random-sampling.html"); 
		mappings.put("Utilities/Intersect lists", "chipster-manual/intersect-lists.html"); 

		mappings.put("aCGH tools/Import from CanGEM", "chipster-manual/import-from-cangem.html");
		mappings.put("aCGH tools/Smooth waves from normalized aCGH data", "chipster-manual/smooth-acgh.html");
		mappings.put("aCGH tools/Call copy number aberrations from aCGH data", "chipster-manual/detect-copy-number-aberrations.html");
		mappings.put("aCGH tools/Plot copy number profiles from called aCGH data", "chipster-manual/plot-cgh-profile.html");
		mappings.put("aCGH tools/Identify common regions from called aCGH data", "chipster-manual/detect-common-copy-number-aberration-regions.html");
		mappings.put("aCGH tools/Cluster called aCGH data", "chipster-manual/cluster-acgh.html");
		mappings.put("aCGH tools/Group tests for called aCGH data", "chipster-manual/stat-acgh.html");
		mappings.put("aCGH tools/Convert called aCGH data from probes to genes", "chipster-manual/convert-cn-probes-to-genes.html");
		mappings.put("aCGH tools/GO enrichment for copy number aberrations", "chipster-manual/pathways-acgh-hyperg-go.html");
		mappings.put("aCGH tools/Match copy number and expression probes", "chipster-manual/match-cn-and-expression-probes.html");
		mappings.put("aCGH tools/Plot profiles of matched copy number and expression", "chipster-manual/plot-cn-induced-expression-profile.html");
		mappings.put("aCGH tools/Test for copy-number-induced expression changes", "chipster-manual/test-for-cn-induced-differential-expression.html");
		mappings.put("aCGH tools/Plot copy-number-induced gene expression", "chipster-manual/plot-cn-induced-gene-expression.html");
		mappings.put("aCGH tools/Fetch probe positions from CanGEM", "chipster-manual/fetch-probe-positions-from-cangem.html");
		mappings.put("aCGH tools/Add cytogenetic bands", "chipster-manual/add-cytobands.html");
		mappings.put("aCGH tools/Count overlapping CNVs", "chipster-manual/count-overlapping-cnvs.html");
		mappings.put("aCGH tools/Update aberration frequencies for called aCGH data", "chipster-manual/calculate-aberration-frequencies.html");
		
		
		mappings.put("Miscellaneous/Multiple sequence alignment", "chipster-manual/seqanal-msa.html");
		mappings.put("Miscellaneous/Phylogenetics", "chipster-manual/seqanal-phylogenetics.html");
	}

	public static String mapToHelppage(OperationDefinition definition) {
		String page = mappings.get(definition.getCategory().getName() + "/" + definition.getDisplayName());
		if (page == null) {
			page = DEFAULT_HELP_PAGE;
		}
		return page;
	}
}
