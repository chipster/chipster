package fi.csc.microarray.client.visualisation.methods;

import java.io.File;
import java.io.IOException;

import javax.swing.JComponent;
import javax.swing.JFrame;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.SimplePDFViewer;

public class PDFViewer extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("application/pdf");
		
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
        SimplePDFViewer viewer = new SimplePDFViewer(true);
        File file = Session.getSession().getDataManager().getLocalFile(data);
        viewer.openFile(file);
		return viewer;
	}

	
	public static void main(String[] args) throws IOException {
        SimplePDFViewer viewer = new SimplePDFViewer(true);
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(viewer);
        frame.pack();
        frame.setVisible(true);
        File file = new File("document.pdf");
        viewer.openFile(file);

	}
}
