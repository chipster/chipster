package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for showing the location of predicted peaks. Peaks cannot overlap. 
 *
 */
public class PeakTrack extends Track {

	private static final int MIN_VISIBLE_PEAK_SIZE = 5;

	private static final int PEAK_SYMBOL_HEIGHT = 10;

	private Collection<RegionContent> peaks = new TreeSet<RegionContent>();
	
	private Color color;

	public PeakTrack(Color color) {

		this.color = color;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (peaks != null) {

			Iterator<RegionContent> iter = peaks.iterator();
			while (iter.hasNext()) {

				RegionContent peak = iter.next();

				if (!getView().requestIntersects(peak.region)) {
					iter.remove();
					continue;
				}

				createDrawable(peak.region.start, peak.region.end, PEAK_SYMBOL_HEIGHT, color, drawables);
			}
		}

		return drawables;
	}

	private void createDrawable(BpCoord startBp, BpCoord endBp, int height, Color c, Collection<Drawable> drawables) {
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;
		
		if (rect.width < MIN_VISIBLE_PEAK_SIZE) {
			rect.width = MIN_VISIBLE_PEAK_SIZE;
		}

		rect.y = getHeight() / 2;
		rect.height = height;

		drawables.add(new RectDrawable(rect, c, c.darker()));
	}

	public void processDataResult(DataResult dataResult) {

		this.peaks.addAll(dataResult.getContents());
	}    
    
    @Override
	public void defineDataTypes() {
		addDataType(DataType.CHROMOSOME);
		addDataType(DataType.START);
		addDataType(DataType.END);
	}
	
	@Override
	public int getHeight() {
		return PEAK_SYMBOL_HEIGHT * 2;
	}
}
