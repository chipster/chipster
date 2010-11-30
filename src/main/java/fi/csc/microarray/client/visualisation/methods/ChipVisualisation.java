package fi.csc.microarray.client.visualisation.methods;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public abstract class ChipVisualisation extends Visualisation {
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
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
		return isTabular(bean) && bean.hasTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);
	}
}
