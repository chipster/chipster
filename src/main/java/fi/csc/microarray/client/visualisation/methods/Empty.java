package fi.csc.microarray.client.visualisation.methods;

import java.util.List;

import javax.swing.JComponent;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;

/**
 * Empty visualisation is always present and spaces reserve space also for the longer names
 * when they aren't present. This is the easiest found way to keep things
 * steady.
 */
public class Empty extends Visualisation {
	
	public Empty(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean datas) {
		return this.getDefaultVisualisation();		
	}
	
	@Override
	public JComponent getVisualisation(List<DataBean> beans) {
		return this.getDefaultVisualisation();		
	}

	@Override
	public boolean canVisualise(DataBean bean) {
		return true;
	}
	
	@Override
	public boolean canVisualise(List<DataBean> beans) {
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
