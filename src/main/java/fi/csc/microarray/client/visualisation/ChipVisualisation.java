package fi.csc.microarray.client.visualisation;

import java.util.LinkedList;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;

public abstract class ChipVisualisation extends Visualisation {
	
	public ChipVisualisation(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public Variable[] getVariablesFor(DataBean dataBean) {
		LinkedList<Variable> vars = new LinkedList<Variable>();
		try {
			extractChips(dataBean, vars);
		} catch (MicroarrayException e) {
			application.reportException(new MicroarrayException("no chips to visualise"));
		}
		return vars.toArray(new Variable[0]);
		
	}

	private void extractChips(DataBean dataBean, LinkedList<Variable> vars) throws MicroarrayException {
		String exprHeader = "/column/";
		String chipHeader = "chip.";
		Table columns = dataBean.queryFeatures("/column/*").asTable();
		
		for (String columnName : columns.getColumnNames()) {
			if (columnName.startsWith(chipHeader)) {
				String chipName = columnName.substring(chipHeader.length());
				String expression = exprHeader + columnName; 
				vars.add(new Variable(chipName, expression));
				// TODO disable log only when data is log-transformed already 
				//vars.add(new Variable("log_2 of " + chipName, "log(" + expression + ")"));
			}
		}
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
