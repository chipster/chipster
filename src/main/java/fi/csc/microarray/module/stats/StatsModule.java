package fi.csc.microarray.module.stats;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.RestrictModifier;
import fi.csc.microarray.databeans.features.stat.LogModifier;
import fi.csc.microarray.databeans.features.stat.NegModifier;
import fi.csc.microarray.databeans.features.table.RowCountProvider;
import fi.csc.microarray.databeans.features.table.TableColumnProvider;
import fi.csc.microarray.module.Module;

public class StatsModule implements Module {

	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("text/tab", false, false, "tab separated values", VisualConstants.ICON_TYPE_TABLE, "tsv");
		manager.plugContentType("text/csv", false, false, "comma separated values", VisualConstants.ICON_TYPE_TABLE, "csv");
	}

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/column", new TableColumnProvider());
		manager.plugFeatureFactory("/rowcount", new RowCountProvider());
	}

	public void plugModifiers(DataManager manager) {
		manager.plugModifier("log", new LogModifier());
		manager.plugModifier("neg", new NegModifier());
		manager.plugModifier("restrict", new RestrictModifier());
	}

}
