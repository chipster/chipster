package fi.csc.microarray.client.visualisation;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;

public abstract class ChipVisualisation extends Visualisation {
	
	public ChipVisualisation(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public Variable[] getVariablesFor(DataBean dataBean) {
		return VisualisationUtilities.getVariablesFiltered(dataBean, "chip.", true);
	}
			
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		boolean isTabular = VisualisationMethod.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
		return isTabular && hasRows(bean) && bean.queryFeatures("/column/chip.*").exists();
	}
	
	protected boolean hasRows(DataBean dataBean) throws MicroarrayException {
		return dataBean.queryFeatures("/rowcount/max/1").asFloat() >= 1;
	}	
}
