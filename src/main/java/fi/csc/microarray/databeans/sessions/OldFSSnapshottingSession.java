package fi.csc.microarray.databeans.sessions;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;


/**

 * @author Aleksi Kallio
 *
 */
public class OldFSSnapshottingSession {

	private static final String METADATA_FILENAME = "snapshot_metadata.txt";

	private static final String ROOT_FOLDER_ID = "0";
	
	private DataManager manager;

	private HashMap<Integer, DataItem> idMap = new HashMap<Integer, DataItem>();
	private HashMap<DataItem, Integer> reversedIdMap = new HashMap<DataItem, Integer>();

	public OldFSSnapshottingSession(DataManager manager) {
		this.manager = manager;
	}
	private File getMetadataFile(File snapshotDir) {
		File metadataFile = new File(snapshotDir.getPath() + File.separator + METADATA_FILENAME);
		return metadataFile;
	}
	public List<DataItem> loadFromSnapshot(File snapshotDirectory, DataFolder parentFolder) throws IOException, MicroarrayException {
		
		File metadataFile = getMetadataFile(snapshotDirectory);
		
		LinkedList<DataItem> newItems = new LinkedList<DataItem>();
		BufferedReader metadataIn = null;
		try {
			
			// load metadata and data
			metadataIn = new BufferedReader(new FileReader(metadataFile));
			for (String line = metadataIn.readLine(); line != null; line = metadataIn.readLine()) {
				
				if (line.startsWith("DATAFOLDER ")) {
					String[] split = line.split(" ");
					String id = split[1];
					DataFolder folder = manager.createFolder("");
					newItems.add(folder);
					mapId(id, folder);
					
//					if (id.equals(ROOT_FOLDER_ID)) {
//						manager.getRootFolder().addChild((FSDataFolder)fetchItem(ROOT_FOLDER_ID));
//					}
					
				} else if (line.startsWith("DATABEAN ")) {
					String[] split = line.split(" ");
					String id = split[1];
					File contentFile = new File(snapshotDirectory.getAbsolutePath() + File.separator + split[2]);
					//DataBean bean = manager.createDataBean("<empty>", contentFile);
					DataBean bean = manager.createDataBean("<empty>", new FileInputStream(contentFile));
					Operation importOperation = new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] {bean});
					bean.setOperation(importOperation); // in case they have no real operation defined later on
					bean.setContentType(manager.guessContentType(contentFile));
					
					newItems.add(bean);
					mapId(id, bean);
					
				} else if (line.startsWith("NAME ")) {
					String[] split = line.split(" ");
					String id = split[1];
					String name = afterNthSpace(line, 2);
					DataItem item = fetchItem(id);
					item.setName(name);
					if (item instanceof DataBean) {
						// update content type now that we have the real filename available
						DataBean bean = (DataBean)item;
						bean.setContentType(manager.guessContentType(name));
					}
					
				} else if (line.startsWith("OPERATION ")) {
					String[] split = line.split(" ");
					String id = split[1];
					DataBean bean = (DataBean)fetchItem(id);
					String[] opName = afterNthSpace(line, 2).split("/");					
					OperationDefinition od = new OperationDefinition(opName[1], new OperationCategory(opName[0]), "", false);
					Operation op = new Operation(od, new DataBean[] { bean });
					bean.setOperation(op);

				} else if (line.startsWith("OPERATION_COLOR ")) {
					String[] split = line.split(" ");
					String id = split[1];
					DataBean bean = (DataBean)fetchItem(id);
					int color = Integer.parseInt(split[2]);
					bean.getOperation().getDefinition().getCategory().setColor(new Color(color));

				} else if (line.startsWith("CHILD ")) {
					String[] split = line.split(" ");
					String childId = split[1];
					String parentId = split[2];
					
					DataFolder parent;
					if (parentId.equals(ROOT_FOLDER_ID)) {
						parent = manager.getRootFolder();
					} else {
						parent = (DataFolder)fetchItem(parentId);
					}
					
					DataItem child = fetchItem(childId);
					parent.addChild(child);
					
				} else if (line.startsWith("NOTES ")) {
					String[] split = line.split(" ");
					String id = split[1];
					String notes = afterNthSpace(line, 2);
					DataBean item = (DataBean)fetchItem(id);
					item.setNotes(notes);

				} else if (line.startsWith("LINK ")) {
					String[] split = line.split(" ");
					Link link= Link.valueOf(split[1]);
					String fromId = split[2];
					String toId = split[3];
					DataBean from = (DataBean)fetchItem(fromId);
					DataBean to = (DataBean)fetchItem(toId);
					from.addLink(link, to);

				} else {
					throw new RuntimeException("metadata error in " + metadataFile.getCanonicalPath() + ": line could not be processed \"" + line + "\"");
				}				
			}
			
		} finally {
			if (metadataIn != null) {
				metadataIn.close();
			}
		}
		
		return newItems;		
	}

	private static String afterNthSpace(String line, int nth) {
		int from = 0;
		for (int i = 0; i < nth; i++) {
			from = line.indexOf(" ", from + 1);
		}
		return line.substring(from + 1);
	}
	
	private void mapId(String id, DataItem item) {
		Integer iid = Integer.parseInt(id);
		idMap.put(iid, item);
		reversedIdMap.put(item, iid);		
	}
	
	private DataItem fetchItem(String id) {
		Integer iid = Integer.parseInt(id);
		return idMap.get(iid);		
	}
}
