package fi.csc.microarray.module.chipster;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.bio.EmbeddedBinaryProvider;
import fi.csc.microarray.databeans.features.bio.IdentifierProvider;
import fi.csc.microarray.databeans.features.bio.NormalisedExpressionProvider;
import fi.csc.microarray.databeans.features.bio.PhenodataProvider;
import fi.csc.microarray.databeans.features.stat.HierarchicalClusterProvider;
import fi.csc.microarray.databeans.features.stat.SomClusterProvider;
import fi.csc.microarray.module.Module;

public class MicroarrayModule implements Module {

	public static final String ANNOTATION_CAT = "Annotation";
	public static final String ANNOTATION_NAME = "Agilent, Affymetrix or Illumina genelist";
	
	public static final String IMPORT_CAT = "Utilities";
	public static final String IMPORT_FROM_ARRAYEXPRESS_NAME = "Import from ArrayExpress";
	public static final String IMPORT_FROM_GEO_NAME = "Import from GEO";
	
	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("application/x-treeview", true, false, "Newick formatted tree from clustering", VisualConstants.ICON_TYPE_TEXT, "tre");
		manager.plugContentType("application/cel", true, false, "Affymetrix CEL", VisualConstants.ICON_TYPE_RAWDATA, "cel");
	}

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/normalised-expression", new NormalisedExpressionProvider());
		manager.plugFeatureFactory("/phenodata", new PhenodataProvider());
		manager.plugFeatureFactory("/identifier", new IdentifierProvider());
		manager.plugFeatureFactory("/embedded-binary-content", new EmbeddedBinaryProvider());
		manager.plugFeatureFactory("/clusters/som", new SomClusterProvider());
		manager.plugFeatureFactory("/clusters/hierarchical", new HierarchicalClusterProvider());
	}

	public void plugModifiers(DataManager manager) {
		// nothing to plug
	}

}
