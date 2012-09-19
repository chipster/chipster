package fi.csc.microarray.client.visualisation.methods;

import java.io.InputStream;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.microarray.client.visualisation.VisualisationFactory;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.exception.MicroarrayException;

public class SamBamViewer extends VisualisationFactory {

	private static int MAX_RECORD_LIMIT = 1024;
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		InputStream in = data.getContentStream(DataNotAvailableHandling.NULL_ON_NA);
		String txt;
		if (in != null) {
			txt = new SamBamUtils().printSamBam(in, MAX_RECORD_LIMIT);
		} else {
			txt = new String(data.getContentBytes(DataNotAvailableHandling.INFOTEXT_ON_NA));
		}
		JTextPane txtPane = TextViewer.makeTxtPane(txt);
		return new JScrollPane(txtPane);
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.getName().endsWith(".bam") || bean.getName().endsWith(".sam"); 
	}
}
