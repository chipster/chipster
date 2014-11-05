package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;

import javax.swing.JComponent;
import javax.swing.JPanel;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class EmptyVisualisation extends Visualisation {
	
	public EmptyVisualisation() {
		// used by VisualisationMethod.getVisualiser()
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		JPanel panel =  new JPanel();
		panel.setBackground(Color.white);
		return panel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return true;
	}
}
