package fi.csc.microarray.client.screen;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;

public class ShowSourceScreen extends ScreenBase implements ClientApplication.SourceCodeListener {

	private JFrame frame = null;
	private JTextArea sourceText = new JTextArea();
	private ClientApplication application = Session.getSession().getApplication();
	private JScrollPane scrollPane;
	
	public JFrame getFrame() {
		// fetch source
		String[] opName = new String[1];
		opName[0] = (String)childScreenParameter;
		try {
			application.fetchSourceFor(opName, this);
		} catch (MicroarrayException e) {
			application.reportException(e);
		}

		// we'll redraw the screen every time
		frame = new JFrame();
		frame.setContentPane(getContentPane());
		frame.setTitle("Source Code");
		return frame;
	}

	private JPanel getContentPane() {
		JPanel contentPane = new JPanel();
		contentPane.setLayout(new BorderLayout());
		sourceText.setText("Please wait while loading source code...");
		sourceText.setEditable(false);
		contentPane.setPreferredSize(new Dimension(800, 400));
		scrollPane = new JScrollPane(sourceText);
		contentPane.add(scrollPane, BorderLayout.CENTER);
		return contentPane;
	}

	public boolean hasFrame() {
		return frame != null;
	}

	public void updateSourceCodeAt(int index, String sourceCode) {
		sourceText.setText(sourceCode);
		sourceText.setCaretPosition(0); // reset scrollpane 
	}

}
