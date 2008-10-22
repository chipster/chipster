package fi.csc.microarray.wizard;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataManager;

public interface WizardContext {

	public DataManager getDataManager();
	public DataSelectionManager getSelectionManager();

	public void importGroup(Collection<ImportItem> datas, String folderName);
	public DataFolder initializeFolderForImport(String string);

	public OperationDefinition locateOperationDefinition(String categoryName, String operationName);
	public void executeOperation(Operation testOp);

	public void visualiseWithBestMethod(FrameType target);
	public void setVisualisationMethod(VisualisationMethod method, List<Variable> variables, List<DataBean> datas, FrameType target );
}
