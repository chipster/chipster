package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.utils.Sequence;

/**
 * Track for showing the reference sequence. Useful only for low zoom levels.
 *
 */
public class SeqTrack extends Track {

	private TreeMap<BpCoord, String> reads = new TreeMap<BpCoord, String>();

	private Color[] charColors = new Color[] { 
			new Color(64, 192, 64, 128), // A
			new Color(64, 64, 192, 128), // C
			new Color(128, 128, 128, 128), // G
			new Color(192, 64, 64, 128), // T
			new Color(192, 192, 192,128) //N
	};

	private long maxBpLength;

	public SeqTrack(View view, DataSource file, long maxBpLength) {

		super(view, file);
		this.maxBpLength = maxBpLength;

	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (reads != null) {

			Iterator<Entry<BpCoord, String>> iter = reads.entrySet().iterator();
			
			while (iter.hasNext()) {

				Entry<BpCoord, String> read = iter.next();

				if (!getView().getBpRegion().contains(read.getKey())) {

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

		if (rect.width > seq.length() * CHAR_WIDTH) {

			drawables.add(new TextDrawable((int) x, rect.y + 9, "" + letter, Color.black));
		}

		Color bg = Color.white;

		if (letter == 'A' || letter == 'a') {
			bg = charColors[0];
		} else if (letter == 'C' || letter == 'c') {
			bg = charColors[1];
		} else if (letter == 'G' || letter == 'g') {
			bg = charColors[2];
		} else if (letter == 'T' || letter == 't') {
			bg = charColors[3];
		} else if (letter == 'N' || letter == 'n') {
			bg = charColors[4];
		}

		nextX = x + increment;
		drawables.add(new RectDrawable(Math.round(x), rect.y - 2,
				Math.round(nextX) - Math.round(x), 10, bg, null));
		x = nextX;

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		
		//Sequence strings has to be cut in pieces to recognise data that we have already
		
		for (RegionContent rc : areaResult.getContents()) {
			
			String seq = ((String) rc.values.get(ColumnType.SEQUENCE));
			
			for ( int i = 0; i < seq.length(); i++ ) {
				
				BpCoord bp = new BpCoord(rc.region.start.bp + i, rc.region.start.chr);
				
				reads.put(bp, "" + seq.charAt(i));
			}
		}
		
		getView().redraw();
	}

	@Override
	public Integer getHeight() {
		if (isVisible()) {
			return 20;
		} else {
			return 0;
		}
	}
	   
    @Override
    public boolean isStretchable() {
        return false;
    }
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() <= maxBpLength);
    }

    @Override
    public Map<DataSource, Set<ColumnType>> requestedData() {
        HashMap<DataSource, Set<ColumnType>> datas = new
        HashMap<DataSource, Set<ColumnType>>();
        datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.SEQUENCE })));
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
	

	@Override
	public String getName() {
		return "Reads";
	}
}
