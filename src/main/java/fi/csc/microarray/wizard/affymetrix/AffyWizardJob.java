package fi.csc.microarray.wizard.affymetrix;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.nexes.wizard.WizardPanelDescriptor;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.table.EditableTable;
import fi.csc.microarray.databeans.features.table.TableBeanEditor;
import fi.csc.microarray.wizard.ResultBlocker;
import fi.csc.microarray.wizard.WizardContext;
import fi.csc.microarray.wizard.WizardParameterBundle;
import fi.csc.microarray.wizard.WizardPlugin;

public class AffyWizardJob {
	
	public LinkedHashMap<String, WizardPanelDescriptor> getPanels() {
		LinkedHashMap<String, WizardPanelDescriptor> panels = new LinkedHashMap<String, WizardPanelDescriptor>();
		panels.put(WizIntroPanelDesc.IDENTIFIER, new WizIntroPanelDesc());
		panels.put(FileGroupPanelDesc.IDENTIFIER, new FileGroupPanelDesc());
		panels.put(NormalisationPanelDesc.IDENTIFIER, new NormalisationPanelDesc());
		panels.put(GroupTestPanelDesc.IDENTIFIER, new GroupTestPanelDesc());
		return panels;
	}
	
	public void execute(WizardContext context, WizardParameterBundle parameters) throws MicroarrayException, IOException {

		// fetch variables
		Map<String, List<String>> fileGroups = parameters.getMappedStrings("fileGroups");
		String test = parameters.getString("test");
		String normalizationMethod = parameters.getString("normalizationMethod");		
		DataManager dm = context.getDataManager();

		// import file groups
		for (String groupName : fileGroups.keySet()) {
			List<String> filegroup = fileGroups.get(groupName);
			List<ImportItem> dataItems = new ArrayList<ImportItem>();
			
			for (Object filename : filegroup) {
				File file = new File(filename.toString());
				ImportItem item = new ImportItem(file);
				item.setType(dm.guessContentType(file.getName()));
				item.setFilename(file.getName());
				
				dataItems.add(item);
			}
						
			context.importGroup(dataItems, groupName);
		}

		// store newly created beans to groups for later use 
		LinkedList<DataBean> datas = new LinkedList<DataBean>();
		HashMap<String, String> groupsForDatas = new HashMap<String, String>();
		for (String group : fileGroups.keySet()) {
			for (DataItem data : dm.getRootFolder().getChildFolder(group).getChildren()) {
				if (data instanceof DataBean) {
					datas.add((DataBean)data);
					groupsForDatas.put(data.getName(), group.substring(WizardPlugin.FILEGROUP_PREFIX.length()));
				}
			}
		}
		
		// create result folder and select it
		context.initializeFolderForImport("results");
		DataFolder resultFolder = dm.getRootFolder().getChildFolder("results");

		// run normalisation
		Operation normOp = new Operation(context.locateOperationDefinition("Normalisation", "Affymetrix"), datas.toArray(new DataBean[0]));
		normOp.setOutputFolder(resultFolder);
		normOp.setParameter("normalization.method", normalizationMethod);	
		ResultBlocker normBlocker = new ResultBlocker(2);	
		normOp.setResultListener(normBlocker);
		context.executeOperation(normOp);
		normBlocker.blockUntilDone();

		// get normalised data
		DataBean normalised = null;
		for (DataBean bean : normBlocker.getResults()) {
			if ("normalized.tsv".equals(bean.getName())) {
				normalised = bean;
				break;
			}
		}
		DataBean phenodata = LinkUtils.retrieveInherited(normalised, Link.ANNOTATION);

		// set groups to phenodata
		// FIXME does this fail when there are more than 9 groups?
		TableBeanEditor tableEditor = new TableBeanEditor(phenodata);
		EditableTable phenodataTable = tableEditor.getEditable();
		for (int i = 0; i < phenodataTable.getRowCount(); i++) {
			String group = groupsForDatas.get(phenodataTable.getValue("original_name", i));
			phenodataTable.setValue("group", i, group);
		}
		tableEditor.write();

		// run statistical test
		String testMethod = test.substring(0, test.indexOf('+'));
		String testAdjustment = test.substring(test.indexOf('+') + 1);

		Operation testOp = new Operation(context.locateOperationDefinition("Statistics", "Several groups tests"), 
				new DataBean[] { normBlocker.getResults().get(0) } );
		testOp.setOutputFolder(resultFolder);
		testOp.setParameter("test", testMethod);	
		testOp.setParameter("p.value.adjustment.method", testAdjustment);
		testOp.setParameter("use.simple.analysis", "yes");
		testOp.setParameter("p.value.threshold", 0.05);
		testOp.setParameter("use.simple.analysis", "yes");
		testOp.setParameter("column", "group");

		ResultBlocker testBlocker = new ResultBlocker(1);	
		testOp.setResultListener(testBlocker);
		context.executeOperation(testOp);
		testBlocker.blockUntilDone();

		// check result
		DataBean testResult = testBlocker.getResults().get(0);
		boolean isEmpty = !testResult.queryFeatures("/column/ ").exists();

		// show empty spreadsheet if that is all we got
		if (isEmpty) {
			context.getSelectionManager().clearAll(true, this);
			context.getSelectionManager().selectSingle(testResult, this);
			context.visualiseWithBestMethod(FrameType.MAIN);

		} else {
			// do hierarchical clustering
			Operation hcOp = new Operation(context.locateOperationDefinition("Clustering", "Hierarchical"), 
					new DataBean[] { testResult } );
			hcOp.setOutputFolder(resultFolder);

			ResultBlocker hcBlocker = new ResultBlocker(2);	
			hcOp.setResultListener(hcBlocker);
			context.executeOperation(hcOp);
			hcBlocker.blockUntilDone();

			// check result				
			DataBean hcResult = hcBlocker.getResults().get(0).getName().contains("hc.txt") ? 
					hcBlocker.getResults().get(0) : hcBlocker.getResults().get(1);

					// visualise HC result
					context.getSelectionManager().clearAll(true, this);
					context.getSelectionManager().selectSingle(hcResult, this);
					context.setVisualisationMethod(VisualisationMethod.HIERARCHICAL, 
							null, 
							context.getSelectionManager().getSelectedDataBeans(), 
							FrameType.MAIN);
		}
	}
}
