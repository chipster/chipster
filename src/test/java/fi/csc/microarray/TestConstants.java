package fi.csc.microarray;


public class TestConstants {

	public static final String TEST_DATA_ROOT = "src/test/resources/";
	
	public static final String FOUR_CHIPS_RESOURCE= TEST_DATA_ROOT + "4_chips_filtered.tsv";
	public static final String FOUR_CHIPS_PHENODATA_RESOURCE = TEST_DATA_ROOT + "4_chips_phenodata.tsv";
	public static final String AFFY_RESOURCE = TEST_DATA_ROOT + "affy_example.cel";
	public static final String HIERARCHICAL_CLUSTERED_AGILENT_RESOURCE = TEST_DATA_ROOT + "agilent-hc.txt";
	public static final String HIERARCHICAL_CLUSTERED_AGILENT_HEATMAP_RESOURCE = TEST_DATA_ROOT + "agilent-hc-normalized.tsv";
	public static final String BIN_AFFY_RESOURCE = TEST_DATA_ROOT + "binary_affy_example.cel";
	public static final String SCATTER_HARDCASE2= TEST_DATA_ROOT + "ei_visualisoidu_2.tsv";
	public static final String SCATTER_HARDCASE1= TEST_DATA_ROOT + "ei_visualisoidu.tsv";
	public static final String HIERARCHICAL_CLUSTERED_RESOURCE = TEST_DATA_ROOT + "hc.tre";
	public static final String HIERARCHICAL_CLUSTERED_HEATMAP_RESOURCE = TEST_DATA_ROOT + "hc-normalized.tsv";
	public static final String HIERARCHICAL_CLUSTERED_ILLUMINA_RESOURCE = TEST_DATA_ROOT + "illumina-hc.tre";
	public static final String HIERARCHICAL_CLUSTERED_ILLUMINA_HEATMAP_RESOURCE = TEST_DATA_ROOT + "illumina-hc-normalized.tsv";
	public static final String CDNA_RESOURCE = TEST_DATA_ROOT + "mouse1.tsv";
	public static final String CLUSTERED_PROFILES_RESOURCE = TEST_DATA_ROOT + "kmeans.tsv";	
	public static final String RESULSET_RESOURCE = TEST_DATA_ROOT + "resultset.tsv";
	public static final String SOM_CLUSTERED_RESOURCE = TEST_DATA_ROOT + "som.tsv";
	
	public static final String DUMMY_ANALYSIS_NAME = "\"Test\"/\"No-op\"";
	public static final String AFFY_NORM_ANALYSIS_NAME = "\"Normalisation\"/\"RMA (Affymetrix)\"";
	public static final String CDNA_NORM_ANALYSIS_NAME = "\"Normalisation\"/\"cDNA\"";
	public static final String AFFY_EXP_ANALYSIS_NAME = "\"Expression\"/\"Expression value filtering\"";
	public static final String LOAD_LIBRARY_ANALYSIS_NAME = "\"Test\"/\"Load\"";
	
	public static final int TIMEOUT_AFTER = 60000;

}
