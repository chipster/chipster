package fi.csc.microarray.client.visualisation.methods;

import java.awt.Font;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class TextViewer extends Visualisation {

	private static long CONTENT_SIZE_LIMIT = 1024*1024*1;
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		byte[] txt = data.getContents(CONTENT_SIZE_LIMIT);

		if (txt != null) {
			JTextPane txtPane = makeTxtPane(new String(txt));
			return new JScrollPane(txtPane);
		}
		return this.getDefaultVisualisation();
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("text/plain", "chemical/x-fasta", "text/wig", "text/bed", "text/fastq", "text/gtf", "text/vcf", "text/qual", "text/mothur-oligos", "text/mothur-names", "text/mothur-groups", "text/sff");
	}
	
	public static JTextPane makeTxtPane(String txt) {
		JTextPane txtPane = new JTextPane();
		txtPane.setFont(Font.decode("Monospaced"));
		txtPane.setText(txt);
		txtPane.setEditable(false);
		return txtPane;
	}
}
