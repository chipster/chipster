package fi.csc.microarray.client.selection;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Selection manager for the rows that are selected from the specific dataset.
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * 
 */
// FIXME should be refactored to not be dataset specific and use IntegratedEntity for all selections
public class IntegratedSelectionManager {

	private ClientApplication client;
	private DataBean data;
	private int[] selectedRows = new int[0];
	private static IntegratedEntity pointSelection; // FIXME remove static and make the whole thing not dataset specific

	public IntegratedSelectionManager(ClientApplication client, DataBean data) {
		this.client = client;
		this.data = data;
	}

	public int[] getSelectionAsRows() {
		return selectedRows.clone();
	}

	public IntegratedEntity getPointSelection() {
		return this.pointSelection;
	}

	public List<String> getSelectionAsIdentifiers() throws MicroarrayException {

		Arrays.sort(selectedRows);

		List<String> names = new ArrayList<String>(selectedRows.length);
		int i = 0;
		for (String name : data.queryFeatures("/identifier").asStrings()) {
			if (Arrays.binarySearch(selectedRows, i) >= 0) {
				names.add(name);
			}
			i++;
		}

		return names;
	}

	public List<String> getSelectedLines() throws Exception {

		List<String> lines = new ArrayList<String>(selectedRows.length + 1);
		BufferedReader original = null;
		original = new BufferedReader(new InputStreamReader(data.getContentStream(DataNotAvailableHandling.EXCEPTION_ON_NA)));
		String line;

		//For binary search
		Arrays.sort(selectedRows);

		for (int i = 0; (line = original.readLine()) != null; i++) {
			if (i == 0 || Arrays.binarySearch(selectedRows, i - 1) >= 0) {
				lines.add(line);
			}
		}

		if (original != null) {
			original.close();
		}

		return lines;
	}
		
	public static DataBean createDataset(Iterable<String> lines, DataBean... sources) throws Exception {
		DataManager dataManager = Session.getSession().getApplication().getDataManager();
		
		DataBean newData = dataManager.createLocalTempDataBean("user_edited.tsv");
		
		// write data
		OutputStream outputStream = dataManager.getContentOutputStreamAndLockDataBean(newData);
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));

		for (String line : lines) {
			writer.write(line);
			writer.newLine();
		}

		writer.flush();
		dataManager.closeContentOutputStreamAndUnlockDataBean(newData, outputStream);
		// TODO get the operation definition from the application
		newData.setOperationRecord(new OperationRecord(new Operation(OperationDefinition.CREATE_DEFINITION, new DataBean[] { newData })));

		
		// set metadata
		DataBean primarySource = sources[0]; // use the first source as the primary source		
		newData.setContentType(primarySource.getContentType());
		for (DataBean data : sources) {
			newData.addLink(Link.MODIFICATION, data);
		}
		dataManager.connectChild(newData, primarySource.getParent());
		
		return newData;
	}

	/**
	 * Normal type of selection.
	 */
	public void setSelection(int[] selection, Object source) {
		selectedRows = selection.clone();
		client.fireClientEvent(new SelectionEvent(data, source));
	}

	/**
	 * Focus type of selection.
	 */
	public void setPointSelection(IntegratedEntity entity, Object source) {
		this.pointSelection = entity;
		client.fireClientEvent(new PointSelectionEvent(data, source));
	}

	public void clearAll(Object source) {
		selectedRows = new int[0];
	}

	public void setSelected(Set<Integer> set, Object source) {

		// To make this little bit more robust, even though null content should
		// be avoided
		set.remove(null);

		int[] indexes = new int[set.size()];

		int i = 0;
		for (Integer index : set) {
			if (index != null) {
				indexes[i++] = index;
			}
		}
		selectedRows = indexes;
		client.fireClientEvent(new SelectionEvent(data, source));
	}
}
