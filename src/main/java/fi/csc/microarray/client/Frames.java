package fi.csc.microarray.client;

import java.awt.Window;
import java.util.LinkedList;

import javax.swing.JFrame;

/**
 * Centralised frame repository for giving frames (JFrame) unified look'n'feel and
 * closing them dangling frames properly at exit.
 * 
 * @author Aleksi Kallio
 *
 */
public class Frames {
	
	public JFrame getMainFrame() {
		return mainFrame;
	}

	private JFrame mainFrame;
	private LinkedList<JFrame> frames = new LinkedList<JFrame>();

	public Frames(JFrame mainFrame) {
		this.mainFrame = mainFrame;
	}
	
	public void registerFrame(JFrame frame) {
		setLocationRelativeToMainFrame(frame);
		frames.add(frame);
	}
	
	public void unregisterFrame(JFrame frame) {
		frames.remove(frame);
	}

	public void setLocationRelativeToMainFrame(Window window) {
		window.setLocationRelativeTo(mainFrame);
	}
}
