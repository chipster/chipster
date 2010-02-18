package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;

public class SeparatorTrack extends Track{

	public SeparatorTrack(View view) {
		super(view, null);
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new LineDrawable(0, 1, getView().getWidth(), 1, Color.gray));
		
		return drawables;				
	}

	public void processAreaResult(AreaResult areaResult) {		
	}
	
	@Override
	public int getMaxHeight(){
		return 3;
	}
	
	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] {}); 
	}
	
	@Override
	public boolean isConcised() {
		return false;
	}
}