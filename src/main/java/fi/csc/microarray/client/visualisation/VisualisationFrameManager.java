package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.Session;

public class VisualisationFrameManager implements PropertyChangeListener{
	
	private VisualisationTaskManager visualisationTaskManager = new VisualisationTaskManager(this);
	
	private InternalVisualisationFrame mainFrame = new InternalVisualisationFrame(FrameType.MAIN);
	private InternalVisualisationFrame sideFrame = new InternalVisualisationFrame(FrameType.SIDE);
	private List<ExternalVisualisationFrame> windows = new ArrayList<ExternalVisualisationFrame>();	
	private JSplitPane split;
	private boolean twinView;
	
	private SimpleInternalFrame frameComponent;

	private VisualisationToolBar toolBar;
	
	public enum FrameType { MAIN, SIDE, WINDOW };
	
	public VisualisationFrameManager(){
		Session.getSession().getApplication().addClientEventListener(this);
		toolBar = new VisualisationToolBar();
	}
	
	public JPanel getFramesPanel(){
		if(frameComponent == null){
			frameComponent = new SimpleInternalFrame("Visualisation");
			frameComponent.setToolBar(toolBar);
			split = new JSplitPane();
			split.setDividerLocation(0.5);
			split.setResizeWeight(0.5);
			this.updateInternalContent();
		} 
		return frameComponent; 
	}
	
	public VisualisationToolBar getVisualisationToolBar(){
		return toolBar;
	}
	
	public Vector<Component> getFocusComponents(){
		Vector<Component> order = new Vector<Component>();
		order.addAll(toolBar.getFocusComponents());		
		return order;
	}
	
	private void updateInternalContent() {
		
		//frameComponent.removeAll();
		int splitLocation = split.getWidth()/2;
		if(twinView){
			splitLocation = split.getDividerLocation();
			split.setLeftComponent(mainFrame.getFrameComponent());
			split.setRightComponent(sideFrame.getFrameComponent());
			frameComponent.setContent(split);			
		} else {
			frameComponent.setContent(mainFrame.getFrameComponent());
		}
		
		//If the program has just started, the splitPane doesn't have width yet
		if(split.getWidth() > 0){
			split.setDividerLocation(splitLocation);
		} else {
			split.setDividerLocation(0.5);	
		}
	}
	
	public void closeAllByType(FrameType type){
		switch (type) {
		case SIDE:
			twinView = false;
			this.updateInternalContent();			
			break;
		case WINDOW:
			for(ExternalVisualisationFrame window: windows){
				window.getFrameComponent().dispose();
			}
			windows.clear();
			break;
		case MAIN:
			throw new IllegalArgumentException("Main visualisation frame can't be closed");
		}	
	}	
	
	public boolean isTwinView(){
		return twinView;
	}

	public void propertyChange(PropertyChangeEvent event) {

		if (event instanceof VisualisationMethodChangedEvent) {
			VisualisationMethodChangedEvent e = (VisualisationMethodChangedEvent) event;
			//logger.debug("VisualisationPanel got VisualisationMethodChangedEvent with method: " + visualisationEvent.getNewMethod());
			
			//Special case: the empty  visualisation is so fast, that showing wait 
			//panel causes only irritating flickering
			if(((VisualisationMethodChangedEvent) event).getNewMethod() != VisualisationMethod.NONE){
				// draw wait panel while executing
				this.showWaitPanel(e.getTarget());
			}
			
			//If visualization is removed (e.g. by opening a new session) in maximized state it becomes difficult to do anything
			if(((VisualisationMethodChangedEvent) event).getNewMethod() == VisualisationMethod.NONE){
				if (toolBar.isMaximised) {
					toolBar.maximiseOrRestoreVisualisation();
				}
			}
			visualisationTaskManager.visualise(e);
		}
	}

	public void showWaitPanel(FrameType frameType) {
		switch(frameType){
		case MAIN:
			mainFrame.showWaitPanel();
			updateInternalContent();
			break;
		case SIDE:
			sideFrame.showWaitPanel();
			twinView = true;
			updateInternalContent();
			break;
		case WINDOW:
			break;
		}
	}

	/**
	 * This method delivers the ready visualisation to the right visualisation frames.
	 * Calls must be done from the EDT because Swing components are manipulated when the 
	 * visualisation is shown.
	 * 
	 * This method SHOULD NOT be called directly outside UpdateGuiRunnable. The 
	 * VisualisationTaskmanager.visualise() prepares visualisations and shows them later 
	 * in right order with this method.
	 * 
	 * @param visualisation
	 * @param target
	 */
	public void showVisualisationComponent(JComponent visualisation, VisualisationMethodChangedEvent e) {
				
		switch (e.getTarget()) {
		case MAIN:
			mainFrame.showVisualisationComponent(visualisation);
			updateInternalContent();
			break;
		case SIDE:
			twinView = true;
			sideFrame.showVisualisationComponent(visualisation);
			updateInternalContent();
			break;
		case WINDOW:						
			e.getTargetFrameInstance().showVisualisationComponent(visualisation);
			break;
		}
	}
	

	public JComponent createVisualisation(VisualisationMethodChangedEvent e) {
		
		switch(e.getTarget()){
		case MAIN:
			return mainFrame.createVisualisation(e);
		case SIDE:
			return sideFrame.createVisualisation(e);			
		case WINDOW:
			ExternalVisualisationFrame window;
			if(e.getTargetFrameInstance() == null || !windows.contains(e.getTargetFrameInstance())){
				window = new ExternalVisualisationFrame();
				e.setTargetFrameInstance(window);
				windows.add(window);
			} else {
				window = (ExternalVisualisationFrame)e.getTargetFrameInstance();
			}
			JComponent visualisation = window.createVisualisation(e);
			return visualisation;
		}
		return null;
	}

	
	public VisualisationFrame getFrame(FrameType target) {
		switch(target){
		case MAIN:
			return mainFrame;
		case SIDE:
			return sideFrame;			
		case WINDOW:						
			//FXIME Support for the others than first
			return windows.get(0);
		}
		return null;
	}
}
