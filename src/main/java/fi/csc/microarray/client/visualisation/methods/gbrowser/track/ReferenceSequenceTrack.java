package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Sequence;

/**
 * Track for showing the reference sequence. Useful only for low zoom levels.
 *
 */
public class ReferenceSequenceTrack extends Track {

	private TreeMap<BpCoord, String> reads = new TreeMap<BpCoord, String>();

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (reads != null) {

			Iterator<Entry<BpCoord, String>> iter = reads.entrySet().iterator();
			
			while (iter.hasNext()) {

				Entry<BpCoord, String> read = iter.next();

				if (!getView().requestContains(read.getKey())) {

					iter.remove();
					continue;
				}
				
				String seq = read.getValue();

				drawables.addAll(getSeqDrawable(read.getKey(), seq, 10));
				drawables.addAll(getSeqDrawable(read.getKey(), Sequence.complement(seq), 0));

			}
		}

		return drawables;
	}

	private Collection<Drawable> getSeqDrawable(BpCoord pos, String seq, int yOffset) {

		Rectangle rect = new Rectangle();

		Collection<Drawable> drawables = getEmptyDrawCollection();

		rect.x = getView().bpToTrack(pos);
		rect.width = getView().bpToTrack(new BpCoord(pos.bp + 1, pos.chr)) - rect.x;

		rect.y = (int) (1 + yOffset);
		rect.height = 10;

		final int CHAR_WIDTH = 7;

		float x = rect.x;
		float increment = (rect.width) / ((float) seq.length());
		float nextX;

		char letter = seq.charAt(0);

		Color bg = Color.white;

		if (letter == 'A' || letter == 'a') {
			bg = GBrowserConstants.charColors[0];
		} else if (letter == 'C' || letter == 'c') {
			bg = GBrowserConstants.charColors[1];
		} else if (letter == 'G' || letter == 'g') {
			bg = GBrowserConstants.charColors[2];
		} else if (letter == 'T' || letter == 't') {
			bg = GBrowserConstants.charColors[3];
		} else if (letter == 'N' || letter == 'n') {
			bg = Color.white;
		}

		nextX = x + increment;
		drawables.add(new RectDrawable(Math.round(x), rect.y - 1,
				Math.round(nextX) - Math.round(x), 10, bg, null));	
		
		if (rect.width > seq.length() * CHAR_WIDTH) {
			
			drawables.add(new TextDrawable((int) x + 1, rect.y + 9, "" + Character.toUpperCase(letter), Color.DARK_GRAY));
		}
		x = nextX;

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {
		
		//Sequence strings has to be cut in pieces to recognise data that we have already
		
		for (Feature rc : dataResult.getFeatures()) {
			
			String seq = ((String) rc.values.get(DataType.SEQUENCE));
			
			if (seq != null) { //when showing negative coordinates
				for ( int i = 0; i < seq.length(); i++ ) {

					BpCoord bp = new BpCoord(rc.region.start.bp + i, rc.region.start.chr);

					reads.put(bp, "" + seq.charAt(i));
				}
			}
		}
	}

	@Override
	public int getTrackHeight() {
		return 20;
	}
    
    @Override
	public void defineDataTypes() {
    	addDataType(DataType.SEQUENCE);
	}

	@Override
	public String getTrackName() {
		return "Reference sequence";
	}
}
