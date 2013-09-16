package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;

/**
 * 3D line that separates real tracks.
 *
 */
public class SeparatorTrack3D extends Track {

	private boolean reversed;

	private static LinkedList<Color> colorSlide = new LinkedList<Color>();
	private static LinkedList<Color> reversedColorSlide = new LinkedList<Color>();

	static {
		Color[] colors = new Color[] {
				new Color(243, 243, 243), 
				new Color(224, 224, 224), 
				new Color(176, 176, 176), 
				new Color(64, 64, 64)
		};
		colorSlide.addAll(Arrays.asList(colors));
		reversedColorSlide.addAll(Arrays.asList(colors));
		Collections.reverse(reversedColorSlide);
	}
	
	
	public SeparatorTrack3D(boolean reversed) {		
		super();
		this.reversed = reversed;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		LinkedList<Color> slide = reversed ? reversedColorSlide : colorSlide;
		
		for (int i = 0; i < slide.size(); i++) {
			drawables.add(new LineDrawable(0, i, getView().getWidth(), i, slide.get(i)));
		}
		
		return drawables;
	}
	
	public void processDataResult(DataResult dataResult) {
		// ignore
	}
	
	@Override
	public int getTrackHeight() {
            return colorSlide.size();
	}

}
