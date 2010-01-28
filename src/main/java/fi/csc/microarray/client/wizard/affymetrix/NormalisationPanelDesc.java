package fi.csc.microarray.client.wizard.affymetrix;

import com.nexes.wizard.*;


public class NormalisationPanelDesc extends WizardPanelDescriptor {
    
    public static final String IDENTIFIER = "NORMALISATION_PANEL";
    
    NormalisationPanel panel;
    
    public NormalisationPanelDesc() {
        panel = new NormalisationPanel();

        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel);
        setButtonPanelType(Wizard.FINISH_NEXT_NAVIGATION_PANEL);
    	//setButtonPanelType(Wizard.NEXT_NAVIGATION_PANEL);
    }
    
    public Object getNextPanelDescriptor() {
    	return GroupTestPanelDesc.IDENTIFIER;
    }
    
    public Object getBackPanelDescriptor() {
        return FileGroupPanelDesc.IDENTIFIER;
    }

    public String getSelected() {
    	return panel.getSelected();
    }
}
