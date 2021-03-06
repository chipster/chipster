package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.session.SessionManager.SessionChangedEvent;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataItemRemovedEvent;

public class VisualisationFrameManager implements PropertyChangeListener, DataChangeListener{
	
	private VisualisationTaskManager visualisationTaskManager = new VisualisationTaskManager(this);
	
	private InternalVisualisationFrame mainFrame = new InternalVisualisationFrame(FrameType.MAIN);
	private InternalVisualisationFrame sideFrame = new InternalVisualisationFrame(FrameType.SIDE);
	private List<ExternalVisualisationFrame> windows = new ArrayList<ExternalVisualisationFrame>();	
	private JSplitPane split;
	private boolean twinView;
	
	private SimpleInternalFrame frameComponent;

	private VisualisationToolBar toolBar;

	private JComponent focusComponent;

	private SwingClientApplication application = (SwingClientApplication) Session.getSession().getApplication();
	
	public enum FrameType { MAIN, SIDE, WINDOW };
	
	public VisualisationFrameManager(){
		toolBar = new VisualisationToolBar();
		application.addClientEventListener(this);
		application.getDataManager().addDataChangeListener(this);
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

		if (focusComponent != null) {
			order.add(focusComponent);
		}
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
			
			// execute visualisation change events
			
			VisualisationMethodChangedEvent e = (VisualisationMethodChangedEvent) event;
			
			//Special case: the empty  visualisation is so fast, that showing wait 
			//panel causes only irritating flickering
			VisualisationMethod method = ((VisualisationMethodChangedEvent) event).getNewMethod();
						
			if(!VisualisationMethod.isDefault(method)){
				// draw wait panel while executing
				this.showWaitPanel(e.getTarget());
			}
			
			//If visualization is removed (e.g. by opening a new session) in maximized state it becomes difficult to do anything
			if (toolBar.isMaximised) {
				toolBar.maximiseOrRestoreVisualisation();
			}
			visualisationTaskManager.visualise(e);
			
		} else if (event instanceof DatasetChoiceEvent) {
			
			// reset visualisation when selection changes
			
			List<DataBean> currentDatas = application .getSelectionManager().getSelectedDataBeans();
			List<DataBean> newDatas = getVisualisedDatas();					
			
			if (currentDatas == null || newDatas == null || // prevent npe 
					// update if selection has changed					
					!(currentDatas.containsAll(newDatas) && newDatas.containsAll(currentDatas))) {
				
				application.setVisualisationMethodToDefault();
			}
		} else if (event instanceof SessionChangedEvent) {					
			if (VisualisationMethod.isDefault(getFrame(FrameType.MAIN).getMethod())) {
				
				// always update default visualisations
				application.setVisualisationMethodToDefault();
			}
		}
	}
	
	private List<DataBean> getVisualisedDatas() {
		return getFrame(FrameType.MAIN).getDatas();
	}

	public void dataChanged(DataChangeEvent e) {
		
		if (e instanceof DataItemRemovedEvent) {
			// reset visualisation if a visualized dataset was removed 
			
			List<DataBean> visualizedDatas = getVisualisedDatas();
			
			if (visualizedDatas != null && visualizedDatas.contains(((DataItemRemovedEvent) e).getDataItem())) {
				application.setVisualisationMethodToDefault();
			}			
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
			this.focusComponent = visualisation;
			application.updateFocusTraversal();
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
			final ExternalVisualisationFrame window;
			if(e.getTargetFrameInstance() == null || !windows.contains(e.getTargetFrameInstance())){
				window = new ExternalVisualisationFrame();
				e.setTargetFrameInstance(window);
				windows.add(window);
			} else {
				window = (ExternalVisualisationFrame)e.getTargetFrameInstance();
			}
			JComponent visualisation = window.createVisualisation(e);
			
			// Other visualisations are removed when the same frame gets the next (often empty) visualisation.
			// In a detached window, we have to listen for window close event for clean up because the frame is
			// disposable and won't be used for other visualisations.
			window.getFrameComponent().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			window.getFrameComponent().addWindowListener(new WindowListener() {
				
				public void windowClosed(WindowEvent e) {
					window.removeVisualisationComponent();
					window.removeVisualiser();
				}
				
				public void windowOpened(WindowEvent e) {}				
				public void windowIconified(WindowEvent e) {}			
				public void windowDeiconified(WindowEvent e) {}			
				public void windowDeactivated(WindowEvent e) {}				
				public void windowClosing(WindowEvent e) {}
				public void windowActivated(WindowEvent e) {}								
			});
			
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
