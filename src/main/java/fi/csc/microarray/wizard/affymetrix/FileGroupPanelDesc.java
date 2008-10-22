package fi.csc.microarray.wizard.affymetrix;

import com.nexes.wizard.*;

import java.util.*;


public class FileGroupPanelDesc extends WizardPanelDescriptor {
    
    public static final String IDENTIFIER = "FILE_GROUPING_PANEL";
    
    FileGroupPanel panel;
    
    public FileGroupPanelDesc() {
        panel = new FileGroupPanel(this);

        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel);
        setButtonPanelType(Wizard.FINISH_NEXT_NAVIGATION_PANEL);
    }
    
    public Object getNextPanelDescriptor() {
    	return NormalisationPanelDesc.IDENTIFIER;
    }
    
    public Object getBackPanelDescriptor() {
        return WizIntroPanelDesc.IDENTIFIER;
    }
    
    
    public void aboutToDisplayPanel() {
    	setNextButton();
    }    

    public Vector getGroups() {
    	return panel.getGroups();
    }
            
    public void setNextButton() {
    	if (panel.isGroupsEmpty()) {
			getWizard().setNextButtonEnabled(false);
			getWizard().setFinishButtonEnabled(false);
		} else {
			getWizard().setNextButtonEnabled(true);
			getWizard().setFinishButtonEnabled(true);
		}
	}
}
