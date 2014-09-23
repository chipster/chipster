package fi.csc.microarray.databeans;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.Strings;

public class HistoryText {
	
	private String[] sourceCodes;
	private DataBean data;
	private ClientApplication application = Session.getSession().getApplication();
	
	public HistoryText(DataBean data) {
		this.data = data;
		this.sourceCodes = null;
	}
	
	/**
	 * Generates a String (usually consisting of many, many rows) about the
	 * history of the selected dataset.
	 * @param param 
	 * 
	 * @return A String to be put in the textarea component of this screen.
	 */
	public String getHistoryText(boolean title, boolean name, boolean date, boolean oper, boolean code, boolean notes, boolean param) {
		if (data == null) {
			return null;
		}
		StringBuffer historyText = new StringBuffer();
		int i = 0;
		for (DataBean listData : MicroarrayModule.getSourcePath(data)) {
			if (title) {
				String titleText = "Step " + (i+1);
				historyText.append(titleText + "\n");
				historyText.append(Strings.repeat("-", titleText.length()) + "\n\n");
			}	
			if (name) {
				historyText.append("Dataset name: " + listData.getName() + "\n");
			}
			if (date) {
				historyText.append("Created " + listData.getDate().toString() + "\n");
			}
			if (oper) {
				OperationRecord operationRecord = listData.getOperationRecord();
				historyText.append("Created with operation: ");
				if (operationRecord != null) {
					historyText.append(operationRecord.getFullName() + "\n");
					if (param) {

						for (ParameterRecord parameterRecord : operationRecord.getParameters()) {
							
							// find out human readable value
							OperationDefinition tool = application .getOperationDefinition(operationRecord.getNameID().getID());
							
							String valueString = null;
												
							if (tool != null) {
								valueString = tool.getHumanReadableParameterValue(parameterRecord);								
							} else {						
								valueString = parameterRecord.getValue();
							}							
							
							historyText.append("Parameter " + parameterRecord.getNameID().getDisplayName() + ": " +
									valueString + "\n");
							
						}
					} else {
                        historyText.append("\n");               
                    }
                    
                    if (code) {
                        String[] sources;
                        try {
                            sources = getSourceCodes();                 
                            if (sources[i] != null) {
                                historyText.append("Operation source code:\n" + Strings.indent(sources[i], 4)  + "\n");
                            } else {
                                historyText.append("Operation source code: <not available>\n");
                            }
                        } catch (MicroarrayException e) {
                            Session.getSession().getApplication().reportException(e);
                        }               
                    }
				} else {
					historyText.append("<last operation unknown>\n");
				}
			}
			if (notes) {
				String notesText = listData.getNotes();
				if (notesText != null) {
					historyText.append("Notes: " + notesText + "\n");
				}
			}
			if (listData != data) {
				// If we aren't yet at the end of the list:
				historyText.append("\n");
			}
			i++;			
		}
		return historyText.toString();
	}
	
	private String[] getSourceCodes() throws MicroarrayException {
		if (sourceCodes == null) {
			
			// make list of wanted source codes
			DataBean[] sourceDatas = MicroarrayModule.getSourcePath(data);
			sourceCodes = new String[sourceDatas.length];
			for (int i = 0; i < sourceDatas.length; i++){
				// might be null, is ok
				sourceCodes[i] = sourceDatas[i].getOperationRecord().getSourceCode();
			}
		}
		return sourceCodes;
	}
}
