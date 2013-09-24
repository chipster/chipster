package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.awt.Color;

/**
 * Data items (like FileLine) have to implement this interface to be able to show up in 
 * ScatterplotTrack. 
 * 
 * @author klemela
 */
public interface ScatterplotValue {
	
	public Float getScatterplotValue();
	public Color getScatterplotColor();
}
