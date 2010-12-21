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
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.CytobandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**  
 * The appearance of a chromosome. Used for high level navigation.
 */
public class CytobandTrack extends Track {

	private static final int THICKNESS = 11;

	private static final int MARGIN = 2;

	private BpCoord maxBp;

	private Collection<RegionContent> bands = new TreeSet<RegionContent>();

	private boolean showText;

	enum Band {

		WHITE("gneg", Color.white), 
		LIGHT_GRAY("gpos25", Color.lightGray),
		MOUSE_LIGHT_GRAY("gpos33", Color.lightGray),
		MID_GRAY("gpos50", Color.gray),
		MOUSE_DARK_GRAY("gpos66", Color.darkGray),
		DARK_GRAY("gpos75", Color.darkGray), 
		BLACK("gpos100", Color.black), 
		GAP("acen", null), 
		OTHER("gvar", Color.white),
		STALK("stalk", null);

		private String id;
		private Color color;

		Band(String id, Color color) {
			this.id = id;
			this.color = color;
		}

		public String getId() {
			return id;
		}

		public Color getColor() {
			return color;
		}
	}

	public Band getBand(String id) {
		for (Band band : Band.values()) {
			if (band.getId().equals(id)) {
				return band;
			}
		}
		return null;
	}

	public CytobandTrack(View view, ChunkDataSource file,
	        Class<? extends AreaRequestHandler> handler, boolean showText) {
		super(view, file, handler);

		this.showText = showText;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (bands != null) {
			

			boolean firstGap = true;
			
			//Iterator instead of for to be able to delete items from the collection
			Iterator<RegionContent> bandRegionIter = bands.iterator();

			while (bandRegionIter.hasNext()) {
								
				RegionContent bandRegion = bandRegionIter.next();
				
				//Remove items that don't belong to this view area
				if (!getView().getBpRegion().intercepts(bandRegion.region)) {
					bandRegionIter.remove();
					continue;
				}

				Band band = getBand((String) bandRegion.values.get(ColumnType.VALUE));
				String text = (String) bandRegion.values.get(ColumnType.ID);

				if (text == null) {
					text = "";
				}

				if (band != null && band.getColor() != null) {
					RectDrawable box = createDrawable(bandRegion.region.start, bandRegion.region.end, band.color);
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

				} else if (band == Band.GAP) {

					int y = (int) MARGIN + 1;

					int sideX = getView().bpToTrack(bandRegion.region.end);
					int cornerX = getView().bpToTrack(bandRegion.region.start);

					if (firstGap) {
						int tmp = sideX;
						sideX = cornerX;
						cornerX = tmp;
						firstGap = false;
					}

					drawables.add(new LineDrawable(sideX, y, cornerX, y + THICKNESS / 2, Color.black));

					drawables.add(new LineDrawable(sideX, y + THICKNESS - 1, cornerX, y + THICKNESS / 2, Color.black));

				} else if (band == Band.STALK) {
					
					Rectangle rect = new Rectangle();

					rect.x = getView().bpToTrack(bandRegion.region.start);
					rect.width = getView().bpToTrack(bandRegion.region.end) - rect.x;

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

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		if (areaResult.content.values.containsKey(ColumnType.METADATA) &&
				((Map<?, ?>)(areaResult.content.values.get(ColumnType.METADATA))).containsKey(CytobandParser.LAST_ROW_OF_CHROMOSOME) && 
				getView().getBpRegion().start.chr.equals(areaResult.content.region.start.chr)) {					

			if (maxBp == null || maxBp.compareTo(areaResult.content.region.end) < 0) {
				maxBp = new BpCoord(areaResult.content.region.end);
			}
		}				

		if (getView().getBpRegion().intercepts(areaResult.content.region)) {

			this.bands.add(areaResult.content);
			getView().redraw();
		}

		//		 this.reads.addAll(result.collection);
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
                ColumnType.ID,
                ColumnType.VALUE })));
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}

	@Override
	public BpCoord getMaxBp(Chromosome chr) {

		return maxBp;
	}
}
