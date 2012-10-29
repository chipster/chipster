package fi.csc.microarray.client.visualisation;

import java.awt.Dimension;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import javax.swing.JComponent;
import javax.swing.JFrame;

import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;

public class ExternalVisualisationFrame extends VisualisationFrame{
	

	private JFrame frameComponent;
	
	public ExternalVisualisationFrame() {
		super(FrameType.WINDOW);
	}
	
	public JFrame getFrameComponent(){
		if(frameComponent == null){
			frameComponent = new JFrame();
			
			frameComponent.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			
			frameComponent.addWindowListener(new WindowListener() {
				
				@Override
				public void windowOpened(WindowEvent arg0) {
				}
				
				@Override
				public void windowIconified(WindowEvent arg0) {
				}
				
				@Override
				public void windowDeiconified(WindowEvent arg0) {
				}
				
				@Override
				public void windowDeactivated(WindowEvent arg0) {
				}
				
				@Override
				public void windowClosing(WindowEvent arg0) {
				}
				
				@Override
				public void windowClosed(WindowEvent e) {
					getVisualisation().removeVisualisation();
				}
				
				@Override
				public void windowActivated(WindowEvent arg0) {
				}
			});
		}
		return frameComponent;
	}
	
	@Override 
	void setTitle(String title){
		getFrameComponent().setTitle(title);
	}

	@Override
	public void setContent(JComponent visualisationComponent) {
		getFrameComponent().add(visualisationComponent);	
		getFrameComponent().setSize(new Dimension(800, 600) );
		getFrameComponent().setVisible(true);
	}

	@Override
	protected void updateContextLinks() {
		// do nothing
	}
}
