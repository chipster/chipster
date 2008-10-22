package fi.csc.microarray.client.visualisation;

import java.awt.BorderLayout;

import javax.swing.JComponent;
import javax.swing.JPanel;

import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;

public class InternalVisualisationFrame extends VisualisationFrame {
	
	private JPanel panel;
	
	public InternalVisualisationFrame(FrameType type){
		super(type);		
	}
	
	public JComponent getFrameComponent(){
		if(panel == null){
			panel = new JPanel(new BorderLayout());
		}
		return panel;
	}

	@Override
	public void setContent(JComponent visualisationComponent) {
		getFrameComponent().add(visualisationComponent, BorderLayout.CENTER);
	}
}
