package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.FlowLayout;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JPanel;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;

/**
 * Empty visualisation is always present and spaces reserve space also for the longer names
 * when they aren't present. This is the easiest found way to keep things
 * steady.
 */
public class Empty extends Visualisation {
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean bean) {
		int selectedDataCount = bean != null ? 1 : 0;
		return generateContextlinkPanel(selectedDataCount);		
	}
	
	private boolean isDataLoaded() {
		return !Session.getSession().getDataManager().databeans().isEmpty();
	}
	
	@Override
	public JComponent getVisualisation(List<DataBean> beans) {
		int selectedDataCount = beans != null ? beans.size() : 0;
		return generateContextlinkPanel(selectedDataCount);		
	}

	@Override
	public boolean canVisualise(DataBean bean) {
		return true;
	}
	
	@Override
	public boolean canVisualise(List<DataBean> beans) {
		return true;
	}
	
	@Override
	public boolean isForSingleData(){
		return true;
	}
	
	@Override
	public boolean isForMultipleDatas(){
		return true;
	}
	
	private JComponent generateContextlinkPanel(int selectedDataCount) {
		
		if (!isDataLoaded()) {
			return super.getDefaultVisualisation(); // show empty screen
		}
		
		JPanel mainPanel = new JPanel();
		mainPanel.setBackground(Color.WHITE);

		JPanel linkPanel = Session.getSession().getPrimaryModule().getContextLinkPanel(selectedDataCount);
		if (linkPanel != null) {
			mainPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
			mainPanel.add(linkPanel);
		}
		
		return mainPanel;
	}

}
