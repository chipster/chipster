/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client.screen;

import javax.swing.JFrame;

/**
 * @author akallio
 *
 */
public interface Screen {
	
	public boolean hasFrame();
	public JFrame getFrame();
	public void setChildScreenParameter(Object parameter);
}
