package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.TreeMap;
import java.util.Map.Entry;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class IntensityTrack extends Track {

	private TreeMap<Region, Float> values = new TreeMap<Region, Float>();
	private long minBpLength;
	private Color color;

	public IntensityTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			ReadInstructions<?> readInstructions, Color c, long maxBpLength) { 

		super(view, file, handler, readInstructions);
		
		this.color = c;
		this.minBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		//Collection<Region> toBeRemoved = new ArrayList<Region>();

		for ( Entry<Region, Float> value : values.entrySet()){

			Region reg = value.getKey();
//
//			if(!reg.intercepts(getView().getBpRegion())){
//				toBeRemoved.add(reg);
//			}

			int x1 = getView().bpToTrack(reg.start);
			int x2 = getView().bpToTrack(reg.end);
			int y2 = (int)getView().getTrackHeight();
			
			int val = (int) Math.min(Math.log(value.getValue() * 10) * 5, getView().getTrackHeight()/4);
			int y1 = (int) (-val + y2);

			drawables.add(new RectDrawable(x1, y1, x2 - x1, y2 - y1, color, null));    			    		    		    		    		
		}

//		for(Region reg : toBeRemoved){
//			values.remove(reg);
//		}

		return drawables;
	}

	@Override
	public void updateData(){

		values.clear();
		super.updateData();
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {
					
		if(areaResult.content.values.containsKey(Content.VALUE)){
			values.put(areaResult.content.region, (Float)areaResult.content.values.get(Content.VALUE));			
			getView().redraw();
		}

//		} else if (areaResult.content instanceof List) {

//			List<List<Object>> reads = (List<List<Object>>) areaResult.content;
//
//			List<Long> ruler = getView().getRulerInfo();
//
//			for (int i = 0; i < ruler.size() - 1; i++){
//
//				long rulerPoint = ruler.get(i);
//
//				Region newRegion = new Region(rulerPoint, ruler.get(i+1));
//				if(!values.containsKey(newRegion)){
//					values.put(newRegion, 0f);
//				}
//			}
//
//			for(List<Object> obj: reads){
//
//				Region reg = new Region(-1, -1);
//				reg.start = (Long)obj.get(areaResult.fileDef.indexOf(Content.BP_START));
//				reg.end = reg.start + ((String)obj.get(areaResult.fileDef.indexOf(Content.SEQUENCE))).length();
//
//				for(Region bin : values.keySet()){
//					if(bin.contains(reg.start)){						
//						values.put(bin, values.get(bin) + reg.getLength() / (float)bin.getLength());
//						break;
//					}
//				}				
//			}			
//			getView().redraw();
//		}
	}
	
	public int getMaxHeight(){
		if(getView().getBpRegion().getLength() > minBpLength){
			return super.getMaxHeight();
		} else {
			return 0;
		}
	}
	
	@Override
	public Collection<Content> getDefaultContents() {
		return Arrays.asList(new Content[] {
				//Content.VALUE
		}); 
	}
	
	@Override
	public boolean isConcised() {
		return true;
	}
}
