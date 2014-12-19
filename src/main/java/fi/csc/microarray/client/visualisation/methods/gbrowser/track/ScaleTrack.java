package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.util.ScaleUtil;

public abstract class ScaleTrack extends Track {
	
	private static final int TOP_MARGIN = 0;
	
	private int height;
	private float minValue;
	private float maxValue;
	private int stepCount;
	private float[] scale;
	
	private boolean xAxisVisible = false;

	public ScaleTrack(int height, float minValue, float maxValue, int stepCount) {
		this.height = height;
		this.minValue = minValue;
		this.stepCount = stepCount;
		setMaxValue(maxValue);				
	}
	
	public ScaleTrack(int height, float minValue, float maxValue) {
		this(height, minValue, maxValue, getDefaultStepCount(height));			
	}
	
	private static int getDefaultStepCount(int height) {
		int stepCount = ScaleUtil.DEFAULT_STEP_COUNT;
		if (height < 40) {
			stepCount = 2;
		}
		return stepCount;
	}

	public ScaleTrack(int height, int stepCount) {
		this(height, 0, 100, stepCount);
	}

	public void setXAxisVisible(boolean visible) {
		xAxisVisible = visible;
	}
	
	public List<Drawable> getScaleDrawables() {
		
		List<Drawable> drawables = new LinkedList<>();
		
		if (getTrackHeight() > 20) {
			
			//Color scaleColor = new Color(220, 220, 220);
			Color scaleColor = Color.black;
			
			final int CHAR_WIDTH = 6;
			final int START_OFFSET = 3;
			final int END_OFFSET = -13;
			int textOffset = START_OFFSET;
			int maxTextWidth = 0;
			
			for (float scaleValue : scale) {
				
				String text = "" + ScaleUtil.format(scaleValue);
				maxTextWidth = Math.max(maxTextWidth, text.length() * CHAR_WIDTH);					 
			}	
			
			int scaleX = super.getVisibleWidth() - maxTextWidth - 4;
			
			for (float scaleValue : scale) {
				
				int lineY = getScaledY(scaleValue);
				int textY = getScaledY(scaleValue) + 6;
				int lineX2 = scaleX - 6;
				String text = "" + ScaleUtil.format(scaleValue);
				
				if (xAxisVisible && scaleValue == 0f) {
					lineX2 = 0;
				}
				
				drawables.add(new LineDrawable(scaleX, lineY, lineX2, lineY, scaleColor));
				drawables.add(new TextDrawable(scaleX + 4, textY + textOffset, text, scaleColor));
				
				textOffset -= (START_OFFSET - END_OFFSET) / scale.length; 
			}	
			
			drawables.add(new LineDrawable(scaleX, getScaledY(minValue), scaleX, getScaledY(maxValue), scaleColor));		
		}
		return drawables;
	}
	
	public int getScaledY(float value) {
		return (int) ((value - minValue) / (maxValue - minValue) * (height - TOP_MARGIN - 2)  + 1);
	}
		
	@Override
	public int getTrackHeight() {
		return height;
	}

	public void setMaxValue(float maxValue) {
		this.maxValue = maxValue;	
		
		scale = ScaleUtil.generateScaleValues(minValue, maxValue, stepCount);
		this.minValue = scale[0];
		this.maxValue = scale[stepCount - 1];
	}
	
	public Integer getMaxY(Collection<Drawable> drawables) {

		int maxY = 0;

		for (Drawable drawable : drawables) {
			maxY = Math.max(drawable.getMaxY() + 1, maxY);            
		}

		return maxY;
	}

	/**
	 * Return current maxValue. Subclasses can override this and calculate a 
	 * new maxValue from the current data.
	 * 
	 * @return
	 */
	public int getMaxTotalCoverage() {
		return (int) maxValue;
	}
	
	public void updateScale() {
		setMaxValue(view.parentPlot.getReadScale().numReads);
	}
	
	@Override
	public void updateLayout() {
		/*
		 * Automatic scale must be calculated for all tracks before drawing any
		 * of them. This is a handy place for such a job.
		 * 
		 * New data may arrive between this calculation and paintComponent(),
		 * which don't fit in this scale, but it will cause a new repaint which
		 * will fix it.
		 */
 
		if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
			view.parentPlot.getReadScale().set(getMaxTotalCoverage());
		}		
	}
}
