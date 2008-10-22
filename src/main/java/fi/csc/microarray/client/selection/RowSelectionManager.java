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

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;

/**
 * Selection manager for the rows that are selected from the specific dataset.
 * 
 * @author klemela
 * 
 */
public class RowSelectionManager {
	/**
	 * Logger for this class
	 */
	// private static final Logger logger =
	// Logger.getLogger(RowSelectionManager.class);
	private ClientApplication client;
	private DataBean data;
	private int[] selectedRows = new int[0];

	public RowSelectionManager(ClientApplication client, DataBean data) {
		this.client = client;
		this.data = data;
	}

	public int[] getSelectedRows() {
		return selectedRows.clone();
	}

	public List<String> getSelectedIdentifiers() throws MicroarrayException {

		Arrays.sort(selectedRows);

		List<String> names = new ArrayList<String>(selectedRows.length);
		int i = 0;
		for (String name : data.queryFeatures("/column/ ").asStrings()) {
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
		original = new BufferedReader(new InputStreamReader(data.getContentByteStream()));
		String line;

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
		DataBean newData = Session.getSession().getApplication().getDataManager().createDataBean("user_edited.tsv");
		
		// write data
		OutputStream outputStream = newData.getContentOutputStreamAndLockDataBean();
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));

		for (String line : lines) {
			writer.write(line);
			writer.newLine();
		}

		writer.flush();
		newData.closeContentOutputStreamAndUnlockDataBean(outputStream);
		newData.setOperation(new Operation(OperationDefinition.USER_MODIFICATION_DEFINITION, new DataBean[] { newData }));

		
		// set metadata
		DataBean primarySource = sources[0]; // use the first source as the primary source		
		newData.setContentType(primarySource.getContentType());
		for (DataBean data : sources) {
			newData.addLink(Link.MODIFICATION, data);
		}
		primarySource.getParent().addChild(newData);
		
		return newData;
	}

	public void setSelected(int[] selected, Object source) {
		selectedRows = selected.clone();
		client.dispatchEvent(new RowChoiceEvent(data, source));
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
		client.dispatchEvent(new RowChoiceEvent(data, source));
	}
}
