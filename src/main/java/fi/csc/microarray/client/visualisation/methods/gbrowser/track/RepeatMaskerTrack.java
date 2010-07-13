package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class RepeatMaskerTrack extends Track{
	
	private long minBpLength;
	private long maxBpLength;
	private Color color;
    
    private Color BACKGROUND = Color.WHITE;
	
	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	
	private Color[] charColors = new Color[] { 
			new Color(64, 192, 64, 128), // A
			new Color(64, 64, 192, 128), // C
			new Color(128, 128, 128, 128), // G
			new Color(192, 64, 64, 128) // T
	};
	
	public RepeatMaskerTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler, 
			Color color, long minBpLength, long maxBpLength){
		
		super(view,file,handler);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		// draw a black rectangle as the background
        drawables.add(new RectDrawable(new Rectangle(0, 0,
                getView().getWidth(), this.getHeight()), BACKGROUND, BACKGROUND));

		if (reads != null) {

			Iterator<RegionContent> iter = reads.iterator();
			while (iter.hasNext()) {

				RegionContent read = iter.next();

				if (!read.region.intercepts(getView().getBpRegion())) {

					iter.remove();
					continue;
				}

				BpCoord startBp = read.region.start;
				BpCoord endBp = read.region.end;

				String seq = ((String) read.values.get(ColumnType.SEQUENCE));
				if (seq != null) {
					seq = seq.trim();
					
				} else {
					seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".substring(0, (int) (endBp.minus(startBp) + 1));
				}
				
				Character ch = seq.charAt(0);
				int count = 0;
				for (int i=0;i<seq.length();i++){
					
					if (((Character.isUpperCase(ch) == true) && (Character.isUpperCase(seq.charAt(i)) == false)) ||
							((Character.isUpperCase(ch) == false) && (Character.isUpperCase(seq.charAt(i)) == true))){
						//add painting
						count = 0;
					}
					else{
						count++;
					}
					ch = seq.charAt(i);
				}

				
				
				//drawables.addAll(getSeqDrawables(startBp, endBp, seq, 1));
				//drawables.addAll(getSeqDrawables(startBp, endBp, complement(seq), 11));
			}
		}
		
		return drawables;
	}

//	private Collection<Drawable> getSeqDrawables(BpCoord startBb, BpCoord endBp, String seq, int yOffset) {
//
//		Rectangle rect = new Rectangle();
//
//		Collection<Drawable> drawables = getEmptyDrawCollection();
//
//		rect.x = getView().bpToTrack(startBb);
//		rect.width = getView().bpToTrack(new BpCoord(endBp.bp + 1, endBp.chr)) - rect.x;
//
//		rect.y = (int) (1 + yOffset);
//		rect.height = 10;
//
//		final int CHAR_WIDTH = 7;
//
//		float x = rect.x;
//		float increment = (rect.width) / ((float) seq.length());
//
//		for (int j = 0; j < seq.length(); j++) {
//
//			char letter = seq.charAt(j);
//
//			if (rect.width > seq.length() * CHAR_WIDTH) {
//
//				drawables.add(new TextDrawable((int) x, rect.y + 8, "" + letter, Color.black));
//			}
//
//			Color bg = Color.white;
//			
//			if (letter == 'A' || letter == 'a') {
//				bg = charColors[0];
//			} else if (letter == 'C' || letter == 'c') {
//				bg = charColors[1];
//			} else if (letter == 'G' || letter == 'g') {
//				bg = charColors[2];
//			} else if (letter == 'T' || letter == 't') {
//				bg = charColors[3];
//			}
//
//			drawables.add(new RectDrawable((int) x, rect.y - 1, (int) increment, 10, bg, null));
//
//			x += increment;
//		}
//
//		return drawables;
//	}
	
	@Override
	public Collection<ColumnType> getDefaultContents() {
		
		return Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE });
	}
	
	@Override
	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		// TODO 
		this.reads.add(areaResult.content);
        getView().redraw();
	}

	@Override
	public boolean isConcised() {
		
		return false;
	}

	@Override
	public boolean isStretchable() {
		
		return false;
	}

	@Override
    public Integer getHeight() {
        if (isVisible()) {
            return 31;
        } else {
            return 0;
        }
    }
	
	@Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength);
    }
	
	
}
