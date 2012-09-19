package fi.csc.microarray.client.visualisation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.tasks.ResultBlocker;
import fi.csc.microarray.client.visualisation.VisualisationFactory.Variable;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.ThreadUtils;

public class VisualisationUtilities {

	private static final ClientApplication application = Session.getSession().getApplication();

	public static DataBean filterBySelection(List<DataBean> datas) {
		try {

			// Doing this with multiple datas isn't pretty, so here is simple solution
			// for single datas

			if (datas.size() == 1) {
				Collection<String> lines = application.getSelectionManager().getSelectionManager(datas.get(0)).getSelectedLines();
				return IntegratedSelectionManager.createDataset(lines, datas.toArray(new DataBean[datas.size()]));
			} else {

				List<String[]> allColumns = new LinkedList<String[]>();

				// Get list of columns names from every dataset
				for (DataBean data : datas) {
					if (application.getSelectionManager().getSelectionManager(data).getSelectionAsRows().length > 0) {

						allColumns.add(data.queryFeatures("/column/*").asTable().getColumnNames());
					}
				}

				List<String> columnOrder = new ArrayList<String>();

				while (allColumns.size() > 0) {

					// Find the longest remaining column name list
					int mostColumns = 0;

					for (String[] columnList : allColumns) {
						if (columnList.length > allColumns.get(mostColumns).length) {
							mostColumns = allColumns.indexOf(columnList);
						}
					}

					// Add columns from the found list to columnOrder if it isn't there already

					for (String col : allColumns.get(mostColumns)) {
						if (!columnOrder.contains(col)) {
							columnOrder.add(col);
						}
					}
					allColumns.remove(mostColumns);
				}

				// Construct new data lines with new column order and from multiple datas
				Map<String, Map<String, String>> values = getSelectedFromMultipleDatas(datas);

				List<String> lines = new ArrayList<String>();

				String newLine = "";

				// Column header
				for (String colName : columnOrder) {
					newLine += colName + "\t";
				}

				if (newLine.endsWith("\t")) {
					newLine = newLine.substring(0, newLine.length() - 1);
				}

				lines.add(newLine);

				// Actual content
				for (String id : values.keySet()) {
					newLine = "";

					Map<String, String> rowValues = values.get(id);

					for (String colName : columnOrder) {
						String value = rowValues.get(colName);
						if (value == null) {
							value = "";
						}
						newLine += value + "\t";
					}

					if (newLine.endsWith("\t")) {
						newLine = newLine.substring(0, newLine.length() - 1);
					}

					lines.add(newLine);
				}

				return IntegratedSelectionManager.createDataset(lines, datas.toArray(new DataBean[datas.size()]));
			}

		} catch (Exception exp) {
			application.reportException(new MicroarrayException("Unable to create user filtered dataset", exp));
			return null;
		}
	}

	public static Map<String, Map<String, String>> getSelectedFromMultipleDatas(List<DataBean> datas) throws Exception {

		// Maps identifiers to row map, where row map maps column names to values
		Map<String, Map<String, String>> lines = new HashMap<String, Map<String, String>>();

		// Collect all rows to row maps, add all columns of duplicate identifiers to same row map

		for (DataBean data : datas) {
			Table columns = data.queryFeatures("/column/*").asTable();

			int[] indexes = application.getSelectionManager().getSelectionManager(data).getSelectionAsRows();

			Arrays.sort(indexes);

			for (int i = 0; columns.nextRow(); i++) {

				if (Arrays.binarySearch(indexes, i) < 0) {
					continue;
				}

				Map<String, String> newColumns = new HashMap<String, String>();

				for (String columnName : columns.getColumnNames()) {
					newColumns.put(columnName, columns.getValue(columnName).toString());
				}

				// TODO should use Feature API for this, but it is not that easy...
				String id = newColumns.containsKey(" ") ? newColumns.get(" ") : newColumns.get("identifier");

				if (!lines.containsKey(id)) {
					lines.put(id, newColumns);
				}

				lines.get(id).putAll(newColumns);
			}
		}

		return lines;
	}

	public static void annotateBySelection(List<DataBean> datas, final String annotationOperationName) {

		try {
			final DataBean filterBySelection = filterBySelection(datas);

			Thread thread = ThreadUtils.getBackgroundThread(new Runnable() {
				public void run() {
					try {

						Operation annotationOperation = new Operation(application.getOperationDefinition(annotationOperationName), new DataBean[] { filterBySelection });
						ResultBlocker opBlocker = new ResultBlocker();
						annotationOperation.setResultListener(opBlocker);
						application.executeOperation(annotationOperation);

					} catch (MicroarrayException e) {
						application.reportException(e);
					}
				}
			});
			thread.start();
		} catch (Exception exp) {
			application.reportException(new MicroarrayException("Unable to collect identifiers for annotation", exp));
		}

	}

	public static Variable[] getVariablesFilteredInclusive(DataBean dataBean, String startsWith, boolean removeStart) {
		String exprHeader = "/column/";

		LinkedList<Variable> vars = new LinkedList<Variable>();
		try {
			Table columns = dataBean.queryFeatures("/column/*").asTable();

			for (String columnName : columns.getColumnNames()) {
				if (columnName.startsWith(startsWith)) {
					String chipName;

					if (removeStart) {
						chipName = columnName.substring(startsWith.length());
					} else {
						chipName = columnName;
					}

					String expression = exprHeader + columnName;
					vars.add(new Variable(chipName, expression));
				}
			}

		} catch (MicroarrayException e) {
			application.reportException(new MicroarrayException("no chips to visualise"));
		}
		return vars.toArray(new Variable[0]);
	}

	public static Variable[] getVariablesFilteredExclusive(DataBean dataBean, Collection<String> columnsToRemove, boolean removeStart) {

		LinkedList<Variable> filteredVars = new LinkedList<Variable>();
		LinkedList<Variable> allVars = new LinkedList<Variable>();
		allVars.addAll(Arrays.asList(getVariablesFilteredInclusive(dataBean, "", false)));
		filteredVars.addAll(allVars);

		String hidden = "chip.";

		for (Variable var : allVars) {
			for (String colToRemove : columnsToRemove) {
				if (var.getName().startsWith(colToRemove)) {
					filteredVars.remove(var);
				}
			}

			if (removeStart && var.getName().startsWith(hidden)) {
				String chipName = var.getName().substring(hidden.length());
				filteredVars.set(filteredVars.indexOf(var), new Variable(chipName, var.getExpression()));
			}

			if (filteredVars.contains(var)) {
				if (var.getName().equals(" ")) {
					filteredVars.set(filteredVars.indexOf(var), new Variable("identifier", var.getExpression()));
				}
			}
		}

		return filteredVars.toArray(new Variable[0]);
	}
}
