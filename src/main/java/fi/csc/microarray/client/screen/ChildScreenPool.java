/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client.screen;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeListener;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.KeyStroke;



/**
 * @author akallio, klemela
 * 
 *
 */
public class ChildScreenPool {

	private HashMap<String, Screen> childScreens = new HashMap<String, Screen>();
	
	private JFrame mainFrame;
	
	public ChildScreenPool(JFrame main){
		mainFrame = main;
	}
	
	/**
	 * The method getFrame() of given child screen should contain at least one JComponent
	 * so that Escape key can be bound to close the screen.
	 * 
	 * @param name
	 * @param childScreen
	 */
	public void put(String name, Screen childScreen) {
		childScreens.put(name, childScreen);				
				
		for(Component comp : childScreen.getFrame().getComponents()){
			if (comp instanceof JComponent) {
				JComponent jcomp = (JComponent) comp;
				
				jcomp.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
						KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "ChildScreenClose");
				
				jcomp.getActionMap().put(
						"ChildScreenClose", new KeyShortcutAction(childScreen));
				break;
			}
		}
	}
	
	private class KeyShortcutAction implements Action {
		private Screen screen;
		
		public KeyShortcutAction(Screen screen){
			this.screen = screen;
		}		

		public void addPropertyChangeListener(PropertyChangeListener listener) {
		}
		
		public Object getValue(String key) {
			return null;
		}
		
		public boolean isEnabled() {
			return true;
		}
		
		public void putValue(String key, Object value) {
		}
		
		public void removePropertyChangeListener(PropertyChangeListener listener) {
		}
		
		public void setEnabled(boolean b) {		
		}
		public void actionPerformed(ActionEvent e) {
			screen.getFrame().dispose();
		}			
	}

	public Screen get(String screenName) {
		return childScreens.get(screenName);
	}
	
	public Iterator<Screen> getScreenIterator(){
		return childScreens.values().iterator();
	}
	
	public void show(String screenName, boolean packed, Object parameter) {
		childScreens.get(screenName).setChildScreenParameter(parameter);
		show(screenName, packed);
	}
	
	public void show(String screenName, boolean packed) {
		JFrame frame = childScreens.get(screenName).getFrame();
		if (packed) {
			frame.pack();
		}
		frame.setLocationRelativeTo(mainFrame);
		frame.setVisible(true);		
	}
	
	public void disposeAll() {
		for (Screen screen : childScreens.values()) {
			if (screen.hasFrame()) {
				screen.getFrame().dispose();
			}
		}
	}

	public boolean hasScreen(String screenName) {
		return childScreens.get(screenName) != null;
	}
}

