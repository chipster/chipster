package fi.csc.microarray.client.dataview;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JScrollPane;
import javax.swing.JViewport;

import org.jgraph.plaf.basic.BasicGraphUI;

/**
 * This class disables autoscroll feature in case it is disables 
 * on the graph's parent scrollbar. The autoscroll feature is great 
 * but it can be super annoying if you don't want it! There was not 
 * any ways to disable autoscroll feature so I had to override the 
 * GraphUI class. Really annoying...
 * 
 * @author mkoski
 *
 */
public class SmartAutoscrollGraphUI extends BasicGraphUI {

	/**
	 * Overrides the helper class which takes care of mouse 
	 * event handling in BasicGraphUI class.
	 * 
	 * @author mkoski
	 *
	 */
	class SmartMouseHandler extends BasicGraphUI.MouseHandler{
		
		/**
		 * Does the same as the mouseDragged method in the super class 
		 * expect this one asks if the autoscroll is disabled.
		 */
		@Override
		public void mouseDragged(MouseEvent e) {
			
			// The following code does autoscroll only if it is 
			// enabled by the parent scrollpane
			
			// Why the heck this need .getParent().getParent and in the autoscroll 
			// method the plain .getParent() is correct?? Don't know but it works...
			if(graph.getParent().getParent() instanceof JViewport 
					&& graph.getParent().getParent().getParent() instanceof JScrollPane){
				JScrollPane scroller = (JScrollPane)graph.getParent().getParent().getParent();
				if(scroller.getAutoscrolls()){
					autoscroll(graph, e.getPoint());
				}
			}
			
			// ** The following code is copied straight from super class
			if (graph.isEnabled()) {
				if (handler != null && handler == marquee)
					marquee.mouseDragged(e);
				else if (handler == null && !isEditing(graph) && focus != null) {
					if (!graph.isCellSelected(focus.getCell())) {
						selectCellForEvent(focus.getCell(), e);
						cell = null;
					}
					if (handle != null)
						handle.mousePressed(e);
					handler = handle;
				}
				if (handle != null && handler == handle)
					handle.mouseDragged(e);
			}
			// ** Super class code ends
		}
	}
	
	/**
	 * Overrides the method which sets mouse listener to 
	 * BasicGraphUI
	 */
	@Override
	protected MouseListener createMouseListener() {
		return new SmartMouseHandler();
	}
	
}
