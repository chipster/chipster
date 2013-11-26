package fi.csc.microarray.client.visualisation.methods;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.BrowsableHtmlPanel;

public class HtmlViewer extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		byte[] html = Session.getSession().getDataManager().getContentBytes(data, DataNotAvailableHandling.EMPTY_ON_NA);
		if (html != null) {
			JTextPane htmlPane = BrowsableHtmlPanel.createHtmlPanel();			
			htmlPane.setText(new String(html));
			return new JScrollPane(htmlPane);
		}
		return this.getDefaultVisualisation();
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("text/html");
	}

}
