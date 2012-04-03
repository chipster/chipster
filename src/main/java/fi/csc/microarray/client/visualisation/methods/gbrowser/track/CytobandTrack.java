package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Cytoband;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**  
 * The appearance of a chromosome. Used for high level navigation.
 */
public class CytobandTrack extends Track {

	private static final int THICKNESS = 11;

	private static final int MARGIN = 2;

	private SortedSet<Cytoband> cbands = new TreeSet<Cytoband>();

	private boolean showText;
	
	final Map<Cytoband.Stain, Color> stainColors = new HashMap<Cytoband.Stain, Color>();
	
	{
		stainColors.put(Cytoband.Stain.GNEG, Color.white);
		stainColors.put(Cytoband.Stain.GPOS25, Color.lightGray);
		stainColors.put(Cytoband.Stain.GPOS33, Color.lightGray);
		stainColors.put(Cytoband.Stain.GPOS50, Color.gray);
		stainColors.put(Cytoband.Stain.GPOS66, Color.darkGray);
		stainColors.put(Cytoband.Stain.GPOS75, Color.darkGray);
		stainColors.put(Cytoband.Stain.GPOS100, Color.black);
		stainColors.put(Cytoband.Stain.GPOS, Color.black);
		stainColors.put(Cytoband.Stain.ACEN, null);
		stainColors.put(Cytoband.Stain.GVAR, Color.lightGray);
		stainColors.put(Cytoband.Stain.STALK, null);
		stainColors.put(Cytoband.Stain.TIP, Color.black);
		stainColors.put(Cytoband.Stain.UNRECOGNIZED, null);
	}

	public CytobandTrack(View view, CytobandDataSource data, boolean showText) {
		super(view, data); 

		this.showText = showText;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (cbands != null) {
			

			boolean firstGap = true;
			
			//Iterator instead of 'for' to be able to delete items from the collection
			Iterator<Cytoband> cbandIter = cbands.iterator();

			while (cbandIter.hasNext()) {
								
				Cytoband cband = cbandIter.next();
				
				//Remove items that don't belong to this view area
				if (!getView().getBpRegion().intersects(cband.getRegion())) {
					cbandIter.remove();
					continue;
				}

				Cytoband.Stain stain = cband.getStain();
				Color stainColor = stainColors.get(stain);
				String text = (String) cband.getBand();

				if (text == null) {
					text = "";
				}

				if (stain != null && stainColor != null) {
					RectDrawable box = createDrawable(cband.getRegion().start, cband.getRegion().end, 
							stainColor);
					drawables.add(box);

					if (showText) {

						final int CHAR_WIDTH = 7;

						int textSpace = 0;
						if (box.x >= 0) {
							textSpace = box.width;
						} else {
							textSpace = box.width + box.x;
						}

						if (textSpace > text.length() * CHAR_WIDTH) {
							int textX = Math.max(0, box.x);
							drawables.add(new TextDrawable(textX, box.y - 2, text, Color.black));
						}
					}
					firstGap = true;

				} else if (stain == Cytoband.Stain.ACEN) {

					int y = (int) MARGIN + 1;

					int sideX = getView().bpToTrack(cband.getRegion().end);
					int cornerX = getView().bpToTrack(cband.getRegion().start);

					if (firstGap) {
						int tmp = sideX;
						sideX = cornerX;
						cornerX = tmp;
						firstGap = false;
					}

					drawables.add(new LineDrawable(sideX, y, cornerX, y + THICKNESS / 2, Color.black));

					drawables.add(new LineDrawable(sideX, y + THICKNESS - 1, cornerX, y + THICKNESS / 2, Color.black));

				} else if (stain == Cytoband.Stain.STALK) {
					
					Rectangle rect = new Rectangle();

					rect.x = getView().bpToTrack(cband.getRegion().start);
					rect.width = getView().bpToTrack(cband.getRegion().end) - rect.x;

					rect.y = (int) MARGIN + THICKNESS / 4;
					rect.height = THICKNESS - THICKNESS / 2;

					drawables.add(new RectDrawable(rect, Color.gray, Color.gray));
				}
			}
		}
		return drawables;
	}

	private RectDrawable createDrawable(BpCoord startBp, BpCoord endBp, Color c) {

		Rectangle rect = new Rectangle();

		// Compensate for border
		rect.x = getView().bpToTrack(startBp) - 1;
		rect.width = getView().bpToTrack(endBp) - rect.x;

		rect.y = (int) MARGIN;
		rect.height = THICKNESS;

		return new RectDrawable(rect, c, Color.black);
	}

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent content : areaResult.getContents()) {
//			if (content.values.containsKey(ColumnType.METADATA) &&
//					((Map<?, ?>)(content.values.get(ColumnType.METADATA))).containsKey(CytobandParser.LAST_ROW_OF_CHROMOSOME) && 
//					getView().getBpRegion().start.chr.equals(content.region.start.chr)) {					
//
//				if (maxBp == null || maxBp.compareTo(content.region.end) < 0) {
//					maxBp = new BpCoord(content.region.end);
//				}
//			}				
//
			
			if (getView().getBpRegion().intersects(content.region)) {
				
				Cytoband cband = (Cytoband) content.values.get(ColumnType.VALUE);
				
				cbands.add(cband);
				
			}						
		}
		getView().redraw();
	}

	@Override
	public Integer getHeight() {
		return showText ? 40 : 20;
	}
	
    @Override
    public boolean isStretchable() {
        return false;
    }

    @Override
    public Map<DataSource, Set<ColumnType>> requestedData() {
    	    	
        HashMap<DataSource, Set<ColumnType>> datas = new
        HashMap<DataSource, Set<ColumnType>>();
        datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.VALUE })));
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
}
