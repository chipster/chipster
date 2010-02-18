package fi.csc.microarray.client.visualisation.methods;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.BrowsableHtmlPanel;

public class HtmlViewer extends Visualisation {

	public HtmlViewer(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		byte[] html = data.getContents();
		if (html != null) {
			JTextPane htmlPane = BrowsableHtmlPanel.createHtmlPanel();			
			String htmlStr = new String(html);
			htmlStr = htmlStr.replace("<head>", "<head>" + VisualConstants.HTML_TABLE_STYLE_INTERNAL);
			htmlStr = htmlStr.replace(" border=\"2\"", "");
			htmlPane.setText(htmlStr);
			return new JScrollPane(htmlPane);
		}
		return this.getDefaultVisualisation();
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("text/html");
	}

}
