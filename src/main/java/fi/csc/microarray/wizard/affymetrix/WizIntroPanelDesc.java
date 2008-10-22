package fi.csc.microarray.wizard.affymetrix;

import com.nexes.wizard.*;


public class WizIntroPanelDesc extends WizardPanelDescriptor {
    
    public static final String IDENTIFIER = "INTRODUCTION_PANEL";
    
    WizIntroPanel panel;
    
    public WizIntroPanelDesc() {
        panel = new WizIntroPanel();

        setPanelDescriptorIdentifier(IDENTIFIER);
        setPanelComponent(panel);
    	setButtonPanelType(Wizard.START_NAVIGATION_PANEL);
    }
    
    public Object getNextPanelDescriptor() {
    	return FileGroupPanelDesc.IDENTIFIER;
    }
    
    public Object getBackPanelDescriptor() {
        return null;
    }
}
