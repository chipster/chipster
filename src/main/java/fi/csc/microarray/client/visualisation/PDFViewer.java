package fi.csc.microarray.client.visualisation;

import javax.swing.JComponent;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class PDFViewer extends Visualisation {

	public PDFViewer(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

}
