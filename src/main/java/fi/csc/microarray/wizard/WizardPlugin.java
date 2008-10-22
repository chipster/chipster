package fi.csc.microarray.wizard;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Vector;

import javax.swing.JOptionPane;

import com.nexes.wizard.Wizard;
import com.nexes.wizard.WizardPanelDescriptor;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.util.ThreadUtils;
import fi.csc.microarray.wizard.affymetrix.AffyWizardJob;
import fi.csc.microarray.wizard.affymetrix.FileGroup;
import fi.csc.microarray.wizard.affymetrix.FileGroupPanelDesc;
import fi.csc.microarray.wizard.affymetrix.GroupTestPanelDesc;
import fi.csc.microarray.wizard.affymetrix.NormalisationPanelDesc;

public class WizardPlugin {
	
	public static String FILEGROUP_PREFIX = "filegroup";

	private SwingClientApplication application;
	private WizardContext context;
	private Wizard wizard = null;
	private LinkedHashMap<String, WizardPanelDescriptor> panels;
	private AffyWizardJob job;
	
	public WizardPlugin(SwingClientApplication application) {
		this.application = application;
		this.context = application;
		
		this.wizard = new Wizard(application.getMainFrame());
		SwingClientApplication.setPlastic3DLookAndFeel(wizard.getDialog());
		this.job = new AffyWizardJob();
		
		this.panels = job.getPanels();		
		for (WizardPanelDescriptor panel : panels.values()) {
			wizard.registerWizardPanel(panel.getPanelDescriptorIdentifier(), panel);
		}

		// make firstly inserted panel to be shown first
        wizard.setCurrentPanel(panels.values().iterator().next().getPanelDescriptorIdentifier());
	}
	
	public void show() {
		int ret = wizard.showModalDialog();		
		
		if (ret == Wizard.FINISH_RETURN_CODE) {
			final Map<String, List<String>> fileGroups = new LinkedHashMap<String, List<String>>();
			
			FileGroupPanelDesc groupDesc = (FileGroupPanelDesc)panels.get(FileGroupPanelDesc.IDENTIFIER);			 
			Vector groups = groupDesc.getGroups();
			for (int i = 0; i < groups.size(); ++i) {
			    FileGroup group = (FileGroup)groups.elementAt(i);
			    List<String> groupList = new LinkedList<String>();
			    fileGroups.put(group.getName(), groupList);
			    
			    ArrayList<File> valueList = group.getFiles();
			    ListIterator<File> iter = valueList.listIterator();
			    while (iter.hasNext()) {
			    	File file = (File)iter.next();
			    	groupList.add(file.getPath());
			    }
			}
			
			try {
				Thread thread = ThreadUtils.getBackgroundThread(new Runnable() {
					public void run() {
						try {
							GroupTestPanelDesc testDesc = (GroupTestPanelDesc)panels.get(GroupTestPanelDesc.IDENTIFIER);			 
							NormalisationPanelDesc normalisationDesc = (NormalisationPanelDesc)panels.get(NormalisationPanelDesc.IDENTIFIER);

							WizardParameterBundle parameters = new WizardParameterBundle();
							parameters.add("fileGroups", fileGroups);
							parameters.add("test", testDesc.getSelected());
							parameters.add("normalizationMethod", normalisationDesc.getSelected());		

							job.execute(context, parameters);
						} catch (Exception e) {
							application.reportException(e);
						}
					}
				});
				thread.start();
				
			} catch (Exception e) {				
				application.reportException(e);
			}
			
			JOptionPane.showMessageDialog(application.getMainFrame(),
					"The wizard is running tasks. Please wait until it is finished.");
		}
	}
}
