package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

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
		
		JPanel panel2 = new JPanel();
		panel2.setBackground(Color.WHITE);
		
		// Initialise context link panel
		JPanel panel = new JPanel();
		panel.setBackground(Color.WHITE);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;

		JLabel label;
		JXHyperlink link;
		
		if (selectedDataCount > 0) {
			label = new JLabel(selectedDataCount + " data(s) selected");
			link = new JXHyperlink(new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {
					System.out.println("go1");
					application.setVisualisationMethod(MicroarrayModule.VisualisationMethods.GBROWSER, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
				}
			});
			link.setText("Open genome browser");

		} else {
			label = new JLabel("No data selected");
			link = new JXHyperlink(new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {
					System.out.println("go2");
					application.selectAllItems();
					System.out.println(application.getSelectionManager().getSelectedDataBeans().size());
					application.setVisualisationMethod(MicroarrayModule.VisualisationMethods.GBROWSER, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
				}
			});
			link.setText("Select all and open genome browser");
		}

		c.insets.set(10, 40, 0, 0);
		panel.add(label, c);
		c.gridy++;
		c.insets.set(10, 50, 0, 0);
		panel.add(link, c);

		panel2.setLayout(new FlowLayout(FlowLayout.LEFT));
		panel2.add(panel);
		
		
		return panel2;

	}

}
