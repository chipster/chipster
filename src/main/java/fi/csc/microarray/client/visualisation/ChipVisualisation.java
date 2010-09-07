package fi.csc.microarray.client.visualisation;

import java.util.Arrays;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule;

public abstract class ChipVisualisation extends Visualisation {
	
	public ChipVisualisation(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public Variable[] getVariablesFor(DataBean dataBean) {
		return VisualisationUtilities.getVariablesFilteredInclusive(dataBean, "chip.", true);
	}
	
	public Variable[] getVariablesMore(DataBean dataBean) {
		
		String[] banList = {
			" ", "symbol", "description", "Probe", "Symbol", "Description", "Chromosome",
			"GenBank", "Cytoband", "UniGene", "PubMed", "GeneOntology", "Pathway", "flag.", "Gene", 
			"Gene.Ontology"
		};
		
		return VisualisationUtilities.getVariablesFilteredExclusive(
				dataBean, Arrays.asList(banList), true);
	}
			
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		boolean isTabular = BasicModule.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
		return isTabular && hasRows(bean) && bean.queryFeatures("/column/chip.*").exists();
	}
	
	protected boolean hasRows(DataBean dataBean) throws MicroarrayException {
		return dataBean.queryFeatures("/rowcount/max/1").asFloat() >= 1;
	}	
}
