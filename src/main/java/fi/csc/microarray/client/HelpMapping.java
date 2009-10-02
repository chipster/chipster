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
		mappings.put("Preprocessing/Impute missing values", "chipster-manual/impute.html");
		mappings.put("Preprocessing/Remove missing values", "chipster-manual/na-omit.html");
		mappings.put("Quality control/Affymetrix basic", "chipster-manual/qc-affy.html");
		mappings.put("Quality control/Affymetrix - using RLE and NUSE", "chipster-manual/qc-affy-rle-nuse.html");
		mappings.put("Quality control/Agilent 1-color", "chipster-manual/qc-agilent-one-color.html");
		mappings.put("Quality control/Agilent 2-color", "chipster-manual/qc-agilent.html");
		mappings.put("Quality control/cDNA", "chipster-manual/qc-cdna.html");
		mappings.put("Quality control/Illumina", "chipster-manual/qc-illumina.html");
		mappings.put("Clustering/Hierarchical", "chipster-manual/cluster-hierarchical.html");
		mappings.put("Clustering/K-Means", "chipster-manual/cluster-kmeans.html");
		mappings.put("Clustering/KNN classification", "chipster-manual/classification-knn.html");
		mappings.put("Clustering/Quality Threshold (QT)", "chipster-manual/cluster-qt.html");
		mappings.put("Clustering/Self-organizing map (SOM)", "chipster-manual/cluster-som.html");
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
		mappings.put("Statistics/Correlate with phenodata", "chipster-manual/stat-correlate-phenodata.html");
		mappings.put("Statistics/Linear modelling", "chipster-manual/stat-linear-modelling.html");
		mappings.put("Statistics/SAM", "chipster-manual/stat-sam.html");
		mappings.put("Statistics/Adjust P-values", "chipster-manual/stat-adjust-p-values.html");
		mappings.put("Normalisation/Affymetrix exon arrays", "chipster-manual/norm-affy-exon.html");
		mappings.put("Normalisation/Affymetrix", "chipster-manual/norm-affy.html");
		mappings.put("Normalisation/Agilent 1-color", "chipster-manual/norm-agilent-1color.html");
		mappings.put("Normalisation/Agilent 2-color", "chipster-manual/norm-agilent.html");
		mappings.put("Normalisation/cDNA", "chipster-manual/norm-cdna.html");
		mappings.put("Normalisation/Illumina", "chipster-manual/norm-illumina.html");
		mappings.put("Normalisation/Illumina - lumi pipeline", "chipster-manual/norm-illumina-lumi.html"); 
		mappings.put("Normalisation/Random effects", "chipster-manual/norm-lme.html");
		mappings.put("Normalization/Normalize to specific samples", "chipster-manual/norm-specific-samples.html");
		mappings.put("Pathways/Bayesian network", "chipster-manual/pathway-bayesian.html");
		mappings.put("Pathways/Boolean network", "chipster-manual/pathway-boolean-bestim.html");
		mappings.put("Pathways/Hypergeometric test for GO", "chipster-manual/pathways-hypergeometric-go.html");
		mappings.put("Pathways/Hypergeometric test for KEGG or PFAM", "chipster-manual/pathways-hypergeometric-kegg.html");
		mappings.put("Pathways/SAFE test for KEGG pathway enrichment", "chipster-manual/pathways-hypergeometric-safe.html");
		mappings.put("Pathways/Gene set test", "chipster-manual/stat-geneset.html");
		mappings.put("Pathways/Protein interactions from IntAct", "chipster-manual/annotate-intact.html");
		mappings.put("Pathways/Associations to Reactome pathways", "chipster-manual/annotate-reactome.html");
		mappings.put("Pathways/Over-representation analysis using ConsensusPathDB", "chipster-manual/annotate-cpdb.html");
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
		mappings.put("Sequence analysis/MSA", "chipster-manual/seqanal-msa.html");
		mappings.put("Sequence analysis/Phylogenetics", "chipster-manual/seqanal-phylogenetics.html");
		mappings.put("Annotation/Affymetrix or Illumina genelist", "chipster-manual/annotate-genelist2html.html");
		mappings.put("Utilities/Word-based query", "chipster-manual/search-queryword.html");
		mappings.put("Utilities/Export GEO's SOFT format", "chipster-manual/export-soft.html");
		mappings.put("Utilities/Export tab2mage format", "chipster-manual/export-tab2mage.html");
		mappings.put("Utilities/Extract genes from clustering", "chipster-manual/extract-genes-from-clustering.html");
		mappings.put("Utilities/Extract genes using a p-value", "chipster-manual/extract-genes-from-stattest.html");
		mappings.put("Utilities/Extract samples from dataset", "chipster-manual/extract-samples-from-dataset.html");
		mappings.put("Utilities/Filter using a column", "chipster-manual/filter-by-column.html");
		mappings.put("Utilities/Average replicate chips", "chipster-manual/average-replicates.html");
		mappings.put("Utilities/Calculate descriptive statistics", "chipster-manual/calculate-descriptive-statistics.html");
		mappings.put("Utilities/Calculate fold change", "chipster-manual/calculate-fold-change.html");
		mappings.put("Utilities/Generate phenodata", "chipster-manual/generate-phenodata.html");
		mappings.put("Utilities/Import from GEO", "chipster-manual/import-from-geo.html");
		mappings.put("Utilities/Merge tables", "chipster-manual/merge-tables.html");
		mappings.put("Utilities/Search by correlation", "chipster-manual/search-correlation.html");
		mappings.put("Utilities/Search by gene name", "chipster-manual/search-queryword.html");
		mappings.put("Utilities/Merge tables", "chipster-manual/merge-tables.html");
		mappings.put("Utilities/Sort samples", "chipster-manual/sort-samples.html");
		mappings.put("Utilities/Delete columns", "chipster-manual/delete-columns.html");
		mappings.put("Utilities/Combine probes to genes", "chipster-manual/combine-probes-to-genes.html");
	}

	public static String mapToHelppage(OperationDefinition definition) {
		String page = mappings.get(definition.getCategory().getName() + "/" + definition.getName());
		if (page == null) {
			page = DEFAULT_HELP_PAGE;
		}
		return page;
	}
}
