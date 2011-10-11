package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;

/**
 * 3D line that separates real tracks.
 *
 */
public class SeparatorTrack3D extends SeparatorTrack {

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
	
	
	public SeparatorTrack3D(View view, long minBpLength, long maxBpLength, boolean reversed) {
		super(view, minBpLength, maxBpLength);
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
	
	@Override
	public Integer getHeight() {
        if (isVisible()) {
            return colorSlide.size();
        } else {
            return 0;
        }
	}

}
