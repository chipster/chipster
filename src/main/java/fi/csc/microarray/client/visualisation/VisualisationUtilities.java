package fi.csc.microarray.client.visualisation;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.selection.RowSelectionManager;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.ThreadUtils;
import fi.csc.microarray.wizard.ResultBlocker;

public class VisualisationUtilities {

	private static final ClientApplication application = Session.getSession().getApplication();

	public static DataBean filterBySelection(List<DataBean> datas) {
		try {
			
			// use the data with most columns as the source for filtering
			int largestDataIndex = 0;
			int largestDataColCount = -1;
			for (int i = 0; i < datas.size(); i++) {
				int colCount = datas.get(i).queryFeatures("/column/*").asTable().getColumnCount();
				if (colCount > largestDataColCount) {
					largestDataColCount = colCount;
					largestDataIndex = i;
				}
			}
			
			Collection<String> lines = application.getSelectionManager().getRowSelectionManager(datas.get(largestDataIndex)).getSelectedLines();
			return RowSelectionManager.createDataset(lines, datas.toArray(new DataBean[datas.size()]));

		} catch (Exception exp) {
			application.reportException(new MicroarrayException("Unable to create user filtered dataset", exp));
			return null;
		}
	}

	public static void annotateBySelection(List<DataBean> datas) {

		try {
			final DataBean filterBySelection = filterBySelection(datas);

			Thread thread = ThreadUtils.getBackgroundThread(new Runnable() {
				public void run() {
					try {

						// run normalisation
						Operation normOp = new Operation(application.locateOperationDefinition(MicroarrayModule.ANNOTATION_CAT, MicroarrayModule.ANNOTATION_NAME), new DataBean[] { filterBySelection });
						ResultBlocker normBlocker = new ResultBlocker(2);
						normOp.setResultListener(normBlocker);
						application.executeOperation(normOp);

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
}
