package fi.csc.microarray.client.visualisation;

import java.awt.Dimension;

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
