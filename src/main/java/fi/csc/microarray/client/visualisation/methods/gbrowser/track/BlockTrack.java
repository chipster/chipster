package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;

public class BlockTrack extends Track{

	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	
	List<Integer> occupiedSpace = new ArrayList<Integer>();

	private Color color;

	private int RESOLUTION = 512;

	private long maxBpLength;

	public BlockTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			FileParser inputParser, Color color, long minBpLength, long maxBpLength)  {
		
		super(view, file, handler, inputParser);		
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
//		Collection<Region> toBeRemoved = new ArrayList<Region>();
		
		occupiedSpace.clear();

		if(reads != null){

			
			
			Iterator<RegionContent> iter = reads.iterator();
			while(iter.hasNext()){

				RegionContent read = iter.next();
				
				Object valueObj = read.values.get(ColumnType.VALUE);
				
				
				if(!read.region.intercepts(getView().getBpRegion())){
					
					iter.remove();
					continue;
				}
				
				if(valueObj == null){
					drawables.add(createDrawable(read.region.start, read.region.end, Color.red));
				} else {
					
					int limited = Math.min((int)(Math.log((Float) valueObj * 100)*50), 255);					
					limited = Math.max(0, limited);
					
					//System.out.println(limited);
					
					Color c = new Color(color.getRed(), color.getGreen(), color.getBlue(), limited);
					int height = 10; //(int) (getView().getTrackHeight() / 2 * limited / 255f);
					drawables.add(createDrawable(read.region.start, read.region.end, height , c));
				}				
			}
		}
				
		return drawables;
	}

	private Drawable createDrawable(BpCoord startBp, BpCoord endBp, Color c){
		return createDrawable(startBp, endBp, 5, c);
	}

	private Drawable createDrawable(BpCoord startBp, BpCoord endBp, int height, Color c){
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;

		int i = 0;

		while(occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x + 1){
			i++;
		}

		int end = rect.x + rect.width;

		if(occupiedSpace.size() > i){
			occupiedSpace.set(i, end);
		} else {
			occupiedSpace.add(end);
		}

		rect.y = (int)(getView().getTrackHeight() - ((i + 1) * (height + 2)));
		rect.height = height;	

		return new RectDrawable(rect, c, null);
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {		

//		if (areaResult.content instanceof List) {
//			List<List<Object>> reads = (List<List<Object>>) areaResult.content;
//
//			for(List<Object> obj: reads){
//
//				//TODO make separate track for database genes and read and remove this hack 
//				long start = (Long)obj.get(areaResult.fileDef.indexOf(Content.BP_START)); 
//				long end;
//				if(areaResult.fileDef.indexOf(Content.BP_END) >= 0){
//					end = (Long)obj.get(areaResult.fileDef.indexOf(Content.BP_END)); 
//				} else {
//					end = start + ((String)obj.get(areaResult.fileDef.indexOf(Content.SEQUENCE))).length();
//				}
//
//				Region reg = new Region(start, end);
//
//				this.reads.add(new RegionValue<Float>(reg, 100f));
//			}			
//			
//		} else 

		if(areaResult.status.concise == isConcised()){

			this.reads.add(areaResult.content);			

			getView().redraw();
		}

		//this.reads.addAll(result.collection);
	}

	private boolean wasLastConsied = true;

	private long minBpLength;
	
	@Override
	public void updateData(){

		if(wasLastConsied != isConcised()){
			reads.clear();
			wasLastConsied = isConcised();
		}
		super.updateData();
	}
	
	@Override
	public int getMaxHeight(){
		if(getView().getBpRegion().getLength() > minBpLength && 
				getView().getBpRegion().getLength() <= maxBpLength){
			
			return super.getMaxHeight();
		} else {
			return 0;
		}
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] {}); 
	}

	@Override
	public boolean isConcised() {
		return getView().getBpRegion().getLength() > 1*1024*1024;
		//return reads.size() > RESOLUTION;
	}
}
