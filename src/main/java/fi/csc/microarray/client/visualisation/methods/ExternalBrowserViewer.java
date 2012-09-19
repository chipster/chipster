package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.io.File;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.VisualisationFactory;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.BrowserLauncher;

/**
 * Opens datasets (files) in external web browser.
 *  
 * @author Aleksi Kallio
 *
 */
public class ExternalBrowserViewer extends VisualisationFactory {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	/**
	 * @return always true (you can dump anything into external web browser)
	 */
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return true;
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		
		// open in browser
		File file = Session.getSession().getDataManager().getLocalFile(data);
		BrowserLauncher.openURL("file://" + file.getAbsolutePath());
		
		// show message to user
		JPanel panel = new JPanel(new BorderLayout());
		JLabel label = new JLabel(data.getName() + " is opened in an external web browser...");
		label.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(label, BorderLayout.CENTER);

		return panel;
	}

}
