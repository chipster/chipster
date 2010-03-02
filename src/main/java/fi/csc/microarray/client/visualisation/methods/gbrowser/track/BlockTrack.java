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

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnDefinition;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Type;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

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
				
				Object valueObj = read.values.get(ColumnType.DESCRIPTION);
				
				
				if(!read.region.intercepts(getView().getBpRegion())){
					
					iter.remove();
					continue;
				}
								
				drawables.add(createDrawable(read.region.start, read.region.end, 10, color));							
			}
		}
				
		return drawables;
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

		if(areaResult.status.concise == this.isConcised() &&
				areaResult.content.values.get(ColumnType.STRAND) == getStrand()){

			this.reads.add(areaResult.content);			

			getView().redraw();
		}
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
		
		return Arrays.asList(new ColumnType[] {		 
				ColumnType.STRAND, ColumnType.DESCRIPTION, ColumnType.VALUE }); 
	}

	@Override
	public boolean isConcised() {
		return false;
	}
}
