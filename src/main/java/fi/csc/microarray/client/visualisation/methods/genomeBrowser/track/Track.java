package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.io.File;
import java.util.Collection;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public abstract class Track implements AreaResultListener {
	
	private View view;
	private File file;
	private boolean isReversed;

	public Track(View view, File file){
		this.view = view;
		this.file = file;
	}
		
	public Track(View view, File file,
			Class<? extends AreaRequestHandler> handler,
			ReadInstructions<?> readInstructions) {
		this(view, file);
		view.getQueueManager().createQueue(file, handler, readInstructions);
	}

	/**
	 * Should be called after Track object is created, but can't be merged to constructor,
	 * because the coming areaResult event could cause call to track object before it's constructed.
	 */
	public void initializeListener(){
		if(file != null){
			view.getQueueManager().addResultListener(file, this);
		}
	}

	public abstract Collection<Drawable> getDrawables();
	
	protected View getView() {
		return view;
	}	

	public void updateData(){		
		if(file != null && this.getMaxHeight() > 0){
			FsfStatus status = new FsfStatus();
			status.clearQueues = true;		
			status.concise = isConcised(); 			
			
			Collection<Content> defCont = getDefaultContents();
			
			//status.file = file; //Done in QueueManager
			//System.out.println(this);
			view.getQueueManager().addAreaRequest(file, 
					new AreaRequest(view.getBpRegion(), defCont, status), true);
		}	
	}
	
	public abstract Collection<Content> getDefaultContents();
	public abstract boolean isConcised();
	
	public Collection<Drawable> getEmptyDrawCollection(){
		return  new LinkedList<Drawable>();
	}
	
	public int getMaxHeight(){
		return Integer.MAX_VALUE;
	}

	public void setReverseStrand(boolean b) {
		this.isReversed = b;
	}
	
	public boolean isReversed() {
		return isReversed;
	}
	
	public boolean isForForwardStrand(AreaResult<RegionContent> areaResult){
		String strand = ((String)areaResult.content.values.get(Content.STRAND));
		return (strand == null) || strand.trim().equals("+") ||
		strand.trim().equals("F");
	}
}
