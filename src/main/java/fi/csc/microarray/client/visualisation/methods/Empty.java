package fi.csc.microarray.client.visualisation.methods;

import java.util.List;

import javax.swing.JComponent;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.Dataset;

public class Empty extends Visualisation {
	public Empty(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public JComponent getVisualisation(Dataset datas) {
		return this.getDefaultVisualisation();		
	}
	
	@Override
	public JComponent getVisualisation(List<Dataset> beans) {
		return this.getDefaultVisualisation();		
	}

	@Override
	public boolean canVisualise(Dataset bean) {
		return true;
	}
	
	@Override
	public boolean canVisualise(List<Dataset> beans) {
		return true;
	}
	
	@Override
	public boolean isForSingleData(){
		return true;
	}
	
	@Override
	public boolean isForMultipleDatas(){
		return true;
	}

}
