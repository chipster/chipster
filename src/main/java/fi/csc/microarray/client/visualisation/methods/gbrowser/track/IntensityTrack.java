package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;

public class IntensityTrack extends Track {

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private long minBpLength;
	private Color color;

	public IntensityTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			FileParser inputParser, Color c, long maxBpLength) { 

		super(view, file, handler, inputParser);

		this.color = c;
		this.minBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		//Collection<Region> toBeRemoved = new ArrayList<Region>();

		for ( RegionContent regCont : values){			

			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end);
			int y2 = (int)getView().getTrackHeight();
			
			int val = (int) Math.min(Math.log((Float)(regCont.values.get(ColumnType.VALUE)) * 100) * 5, 
					getView().getTrackHeight()/4);
			int y1 = (int) (-val + y2);

			drawables.add(new RectDrawable(x1, y1, x2 - x1, y2 - y1, color, null));
			
		}
		

		return drawables;
	}

	@Override
	public void updateData(){

		values.clear();
		super.updateData();
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		if(areaResult.content.values.get(ColumnType.STRAND) == getStrand() && 
				areaResult.content.values.get(ColumnType.VALUE) != null){
			
			values.add(areaResult.content);			
			getView().redraw();
		}
	}

	public int getMaxHeight(){
		if(getView().getBpRegion().getLength() > minBpLength){
			return super.getMaxHeight();
		} else {
			return 0;
		}
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] {
				//Content.VALUE
		}); 
	}

	@Override
	public boolean isConcised() {		
		return true;
	}
}
