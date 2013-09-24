package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ScatterplotFileLineConversion;

public class BedTrackGroup extends TrackGroup {
	
	private RegionTrack regionTrack;
	private ScatterplotTrack scatterplotTrack;
	private ScatterplotFileLineConversion dataThread;

	public BedTrackGroup(GBrowserView view, ScatterplotFileLineConversion dataThread, String title) {
		super(view);
		this.setSettingsEnabled(true);
		this.dataThread = dataThread;
		
		this.setName(title);

		addTracks();	
	}
	
	private void addRegionTrack() {
		SeparatorTrack separator = new SeparatorTrack(Color.white, 20);
		separator.setView(view);
		addTrack(separator);
		regionTrack = new RegionTrack(GBrowserConstants.BED_COLOR);
		regionTrack.setView(view);
		regionTrack.addDataThread(dataThread);
		addTrack(regionTrack);
	}
	
	private void addScatterplotTrack() {
		Float min = dataThread.getMinScatterplotValue();
		Float max = dataThread.getMaxScatterplotValue();
		
		if (min == null) {
			min = 0f;
		}
		
		if (max == null) {
			max = 1f;
		}
		
		scatterplotTrack = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 100, min, max, DataType.VALUE);
		scatterplotTrack.setMinSize(4);
		scatterplotTrack.setXAxisVisible(true);
		scatterplotTrack.setView(view);
		scatterplotTrack.addDataThread(dataThread);		
		addTrack(scatterplotTrack);
	}
	
	@Override
	public void addTracks() {		
		tracks.clear();

		getStatusAnimation().addDataThread(dataThread);

		if (!isMinimized()) {
			if (isShowMore()) {
				addScatterplotTrack();
			} else {
				addRegionTrack();
			}
		}
	}
	
	@Override
	public boolean isShowMorePossible() {
		return !isMinimized();
	}
	
	@Override
	public String getShowMoreName() {
		return "Show score";
	}
}
