package fi.csc.microarray.client.wizard.affymetrix;

import com.nexes.wizard.*;


public class GroupTestPanelDesc extends WizardPanelDescriptor {
    
    public static final String IDENTIFIER = "GROUPTEST_PANEL";

	GroupTestPanel panel;

	public GroupTestPanelDesc() {
		panel = new GroupTestPanel();

		setPanelDescriptorIdentifier(IDENTIFIER);
		setPanelComponent(panel);
		setButtonPanelType(Wizard.FINISH_NAVIGATION_PANEL);
	}

	public Object getNextPanelDescriptor() {
		return FINISH;
	}

	public Object getBackPanelDescriptor() {
		return NormalisationPanelDesc.IDENTIFIER;
	}

	public String getSelected() {
		return panel.getSelected();
	}

	public void aboutToDisplayPanel() {
	}

}
