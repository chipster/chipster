package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class CytobandTrack extends Track{
	
	private static final int THICKNESS = 10;

	private static final int MARGIN = 2;

	private Collection<RegionContent> bands = new TreeSet<RegionContent>();

	List<Integer> occupiedSpace = new ArrayList<Integer>();

	private boolean showText;
	
	enum Band {
		
		WHITE("gneg", Color.white), 
		LIGHT_GRAY("gpos25", Color.lightGray), 
		MID_GRAY("gpos50", Color.gray), 
		DARK_GRAY("gpos75", Color.darkGray), 
		BLACK("gpos100", Color.black), 
		GAP("acen", null), 
		OTHER("gvar", Color.white);
		
		private String id;
		private Color color;

		Band(String id, Color color) 
		{
			this.id = id;
			this.color = color;
		}
		
		public String getId(){
			return id;
		}
		
		public Color getColor(){
			return color;
		}
	}
	
//	private class BandRegion implements Comparable<BandRegion>{
//		public Band band;
//		public Region region;
//		
//		public BandRegion(Band band, Region region){
//			this.band = band;
//			this.region = region;
//		}
//		
//		public int compareTo(BandRegion o) {
//			return region.compareTo(o.region);
//		}		
//	}
	
	public Band getBand(String id){
		for (Band band : Band.values()){
			if(band.getId().equals(id)){
				return band;
			}
		}
		return null;
	}

	public CytobandTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			FileParser inputParser, boolean showText) {
		super(view, file, handler, inputParser);		
		
		this.showText = showText;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		occupiedSpace.clear();

		if(bands != null){

			boolean firstGap = true;
			
			for (RegionContent	bandRegion : bands){
				
				Band band = getBand((String)bandRegion.values.get(ColumnType.VALUE));
				String text = (String)bandRegion.values.get(ColumnType.ID);
				
				if(text == null){
					text = "";
				}				

				if(band != null && band.getColor() != null){
					RectDrawable box = createDrawable(bandRegion.region.start, bandRegion.region.end, band.color);
					drawables.add(box);
					
					if(showText){
						
						final int CHAR_WIDTH = 7;
						
						int textSpace = 0;
						if(box.x >= 0){
							textSpace = box.width;			
						} else {
							textSpace = box.width + box.x;
						}
						
						if(textSpace > text.length() * CHAR_WIDTH){
							
							int textX = Math.max(0, box.x);
												
							drawables.add(new TextDrawable(textX, box.y - 2, 
								text, Color.black));
							
//							drawables.add(new TextDrawable(textX, box.y + 5, 
//									text, Color.black));
						}
					}
					firstGap = true;
				} else if (band == Band.GAP) {
		
					int y = (int)getMaxHeight() - (THICKNESS + MARGIN);					
											
					int sideX = getView().bpToTrack(bandRegion.region.end);
					int cornerX = getView().bpToTrack(bandRegion.region.start);
					
					if(firstGap){
						int tmp = sideX;
						sideX = cornerX;
						cornerX = tmp;
						firstGap = false;
					} 					
					
					drawables.add(new LineDrawable(
							sideX, y, 
							cornerX, y + THICKNESS / 2 , Color.black));
					
					drawables.add(new LineDrawable(
							sideX, y + THICKNESS, 
							cornerX, y + THICKNESS / 2, Color.black));
				}
			}
		}
		return drawables;
	}


	private RectDrawable createDrawable(BpCoord startBp, BpCoord endBp, Color c){
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;

		rect.y = (int)(getMaxHeight() - (THICKNESS + MARGIN));
		rect.height = THICKNESS;

		return new RectDrawable(rect, c, Color.black);
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

//		if (areaResult.content instanceof List) {
//			List<List<Object>> reads = (List<List<Object>>) areaResult.content;
//
//			for(List<Object> obj: reads){
//				Region reg = new Region(
//						(Long)obj.get(areaResult.fileDef.indexOf(Content.BP_START)),
//						(Long)obj.get(areaResult.fileDef.indexOf(Content.BP_END)));
//				
//
//				this.bands.add(new BandRegion(band, reg));
//			}
//		}
			
			this.bands.add(areaResult.content);
		
		getView().redraw();
		//this.reads.addAll(result.collection);
	}
	
	

	@Override
	public void updateData(){

		//bands.clear();
		super.updateData();
	}
	
	@Override
	public int getMaxHeight(){
		return showText ? 40 : 20;
	}
	
	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] {
				ColumnType.VALUE,
				ColumnType.ID
		}); 
	}
	
	@Override
	public boolean isConcised() {
		return false;
	}
}
