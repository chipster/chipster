package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Class gets the data values from the iterables and scales them to the scales that
 * are also decided here.
 * 
 * @author Petri Klemel�
 */
public class DataModel {
	
	private static final Logger logger = Logger
	.getLogger(DataModel.class);
	
	private Color lineColor = Color.DARK_GRAY;

	/**
	 * Class handles the mapping from the scaled values to colors used in visualisation.
	 * According the Visualisation researches the rainbow gradient isn't very optimal, 
	 * because of the unclear ordering of the colors and human fixation to recognize
	 * gradient as separate colors with names. However, other possibilities aren't tried yet.
	 * 
	 * @author Petri Klemela
	 */
	public class ColorModel {
		
		
		//About 0.68
		final float BLUE_HUE = Color.RGBtoHSB(0, 0, 255,new float[3])[0];
		
		/**
		 * @param value
		 *            between 0 and 1.0
		 * @return corresponding color
		 */
		public Color getColorFor(float value) {
			value = BLUE_HUE - BLUE_HUE * value;
			return new Color(Color.HSBtoRGB(value, 1, 1));
		}

	}
	
	private NumberFormat numberFormat = NumberFormat.getInstance();

	private ColorModel colorModel = new ColorModel();

	private Drawable[] points;

	//How many lines are shown in the grid for each axis
	private static final int LINE_COUNT = 4;

	private float[] cScale;
	
	public DataModel(){
		numberFormat.setMaximumFractionDigits(2);
	}
	
	public ColorModel getColorModel() {
		return colorModel;
	}


	public void setData(Iterable<String> identifiers, Iterable<Float> xValues, Iterable<Float> yValues,
			Iterable<Float> zValues, Iterable<Float> cValues) {
		
		
		
		List<Drawable> points = new ArrayList<Drawable>();

		Iterator<String> identIter = identifiers.iterator();
		Iterator<Float> xIter = xValues.iterator();
		Iterator<Float> yIter = yValues.iterator();
		Iterator<Float> zIter = zValues.iterator();
		Iterator<Float> cIter = cValues.iterator();
		
		long time = 0;
		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();
		//TODO very slow, 25 secs
		float[] minMaxX = this.findMinMaxValues(xIter);
		float[] minMaxY = this.findMinMaxValues(yIter);
		float[] minMaxZ = this.findMinMaxValues(zIter);
		float[] minMaxC = this.findMinMaxValues(cIter);
		
		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();

		float[] xScale = this.generateScaleValues(minMaxX[0], minMaxX[1]);
		float[] yScale = this.generateScaleValues(minMaxY[0], minMaxY[1]);
		float[] zScale = this.generateScaleValues(minMaxZ[0], minMaxZ[1]);
		cScale = this.generateScaleValues(minMaxC[0], minMaxC[1]);

		points.addAll(this.generateScaleLines(xScale, yScale, zScale));

		// reset iterators
		xIter = xValues.iterator();
		yIter = yValues.iterator();
		zIter = zValues.iterator();
		cIter = cValues.iterator();

		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();
		
		//TODO very slow, 25 secs with 22000 rows
		xIter = xValues.iterator();
		yIter = yValues.iterator();
		zIter = zValues.iterator();
		cIter = cValues.iterator();
		for (int i = 0; identIter.hasNext() && xIter.hasNext() && yIter.hasNext() && 
				zIter.hasNext() && cIter.hasNext(); i++) {
			String identifier = identIter.next();
			float x = xIter.next();
			float y = yIter.next();
			float z = zIter.next();
			float c = cIter.next();
			float[] scaled = this.convertToScaled(xScale, yScale, zScale,
					cScale, new float[] { x, y, z, c });

			Color color = colorModel.getColorFor(scaled[3]);

			points.add(new DataPoint(scaled[0], scaled[1], scaled[2], color,
					0.01, identifier, i));
		}
				
		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();
		
		this.points = new Drawable[points.size()];
		for (int i = 0; i < points.size(); i++) {
			this.points[i] = points.get(i);
		}
		
		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();
	}

	private float[] findMinMaxValues(Iterator<Float> iter) {
		float x;
		float min = Float.MAX_VALUE;
		float max = Float.MIN_VALUE;
		while (iter.hasNext()) {
			x = iter.next();
			if (x > max) {
				max = x;
			}
			if (x < min) {
				min = x;
			}
		}

		return new float[] { min, max };
	}

	public Drawable[] getDataArray() {
		return points;
	}

	private float[] convertToScaled(float[] xScale, float[] yScale,
			float[] zScale, float[] cScale, float[] original) {		
		
		return new float[] { 
				(original[0] - xScale[0]) / (xScale[LINE_COUNT - 1] - xScale[0]),
				(original[1] - yScale[0]) / (yScale[LINE_COUNT - 1] - yScale[0]),
				(original[2] - zScale[0]) / (zScale[LINE_COUNT - 1] - zScale[0]),
				(original[3] - cScale[0]) / (cScale[LINE_COUNT - 1] - cScale[0]) };
	}
	
	/**
	 * Generates a String triplet divided in to table from the coordinates to be used with
	 * Line-class.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	private String[] coordToStringTable(float x, float y, float z){
		return new String[]{"(", "" + x, ", ","" + y, ", ", "" + z, ")" };
	}

	private List<Drawable> generateScaleLines(float[] xValues, float[] yValues,
			float[] zValues) {
		
		//Colors for coordinate triplet: (xx.xx, yy.yy, zz,zz)
		final Color[] tripletColors =
			new Color[] { Color.WHITE, Color.RED, Color.WHITE,
				Color.GREEN, Color.WHITE, Color.BLUE, Color.WHITE };	

		//Coordinate axles
		List<Drawable> lines = new ArrayList<Drawable>();
		lines.add(new Line(0, 0, 0, 1, 0, 0, Color.RED, 
				coordToStringTable(xValues[LINE_COUNT-1], yValues[0], zValues[0]), tripletColors));				
		lines.add(new Line(0, 0, 0, 0, 1, 0, Color.GREEN, 
				coordToStringTable(xValues[0], yValues[LINE_COUNT-1], zValues[0]), tripletColors));
		lines.add(new Line(0, 0, 0, 0, 0, 1, Color.BLUE, 
				coordToStringTable(xValues[0], yValues[0], zValues[LINE_COUNT-1]), tripletColors));
		
		//And their end arrows
		lines.add(new Line(1,0,0, 0.99, 0.005, -0.01, Color.RED, ""));
		lines.add(new Line(1,0,0, 0.99, -0.01, 0.005, Color.RED, ""));
		
		lines.add(new Line(0,1,0, 0.005, 0.99, -0.01, Color.GREEN, ""));
		lines.add(new Line(0,1,0, -0.01, 0.99, 0.005, Color.GREEN, ""));
		
		lines.add(new Line(0,0,1, -0.01, 0.005, 0.99, Color.BLUE, ""));
		lines.add(new Line(0,0,1, 0.005, -0.01, 0.99, Color.BLUE, ""));
		

		//To enable modification of texts to show all coordinates in the corners
		String[] xText;
		String[] yText;
		String[] zText;
		Color[] xTextColor = new Color[] { Color.WHITE };
		Color[] yTextColor = new Color[] { Color.WHITE };
		Color[] zTextColor = new Color[] { Color.WHITE };

		for (int i = 1; i < LINE_COUNT; i++) {
			
			float pos = (float)i / (LINE_COUNT - 1); 

			if (i != LINE_COUNT - 1) {
				//Normal lines show only one coordinate
				xText = new String[] { numberFormat.format(xValues[i])};
				yText = new String[] { numberFormat.format(yValues[i])};
				zText = new String[] { numberFormat.format(zValues[i])};
				

				xTextColor = new Color[] { Color.RED };
				yTextColor = new Color[] { Color.GREEN };
				zTextColor = new Color[] { Color.BLUE };
			} else {
				//Mark texts disabled ( the last lines )
				xText = null;
				yText = null;
				zText = null;
				
				//Colors for coordinate triplet: (xx.xx, yy.yy, zz,zz)
				xTextColor = tripletColors; 
				yTextColor = tripletColors;
				zTextColor = tripletColors;				
			}
			//To understand following  method calls it's not a bad idea to draw
			//a unit cube and mark coordinates of the corners to it
			/*
			 * 				Y
			 * 	(0, 1, 1) *	|\
			 * 				| \                       
			 * 				|  \
			 * 				|	\___________________ * (1, 1, 0) 
			 * 				|   | * (0, 1, 0)		| 
			 * 				|   |					|
			 *  			|   |					|
			 * 				|   | 					|
			 * 				|   | 					|
			 * 				|   | 					|
			 * 				|   | 					|
			 * 				|   |  					|
			 * 				|   |__________________>|X * (1, 0, 0)						 
			 *				|  / * (0, 0, 0)		\
			 * 				| /						 \
			 * 	(0, 0, 1) *	|/________________________\ * (1, 0, 1)
			 * 				Z
			 */

			// X-scale 
			
			lines.add(new Line(pos, 0, 0, pos, 1, 0, lineColor, xText, xTextColor));
			if (xText == null) { // All cordinates for the corner
				xText = coordToStringTable(xValues[i], yValues[0], zValues[i]);
			}
			lines.add(new Line(pos, 0, 0, pos, 0, 1, lineColor, xText, xTextColor));

			// Y-scale 
			
			lines.add(new Line(0, pos, 0, 0, pos, 1, lineColor, yText, yTextColor));
			if (yText == null) {// All cordinates for the corner
				yText = coordToStringTable(xValues[i], yValues[i], zValues[0]);
			}
			lines.add(new Line(0, pos, 0, 1, pos, 0, lineColor, yText, yTextColor));

			// Z-Scale 
			
			lines.add(new Line(0, 0, pos, 1, 0, pos, lineColor, zText, zTextColor));
			if (zText == null) {// All cordinates for the corner
				zText = coordToStringTable(xValues[0], yValues[i], zValues[i]);
			}
			lines.add(new Line(0, 0, pos, 0, 1, pos, lineColor, zText, zTextColor));
		}

		return lines;

	}

	private float[] generateScaleValues(float min, float max) {

		float[] scaleMinMax = this.getLineDistance(min, max);

		float[] values = new float[LINE_COUNT];
		for (int i = 0; i < LINE_COUNT; i++) {
			values[i] = scaleMinMax[0] + (scaleMinMax[1] - scaleMinMax[0])
					/ (LINE_COUNT - 1) * i;
		}

		return values;
	}

	private float[] getLineDistance(float minValue, float maxValue) {

		float preferredScale = Math.abs(maxValue - minValue) / (LINE_COUNT - 1);
		float decimalFactor = this.calculateDecimalFactor(preferredScale);
		float scaleMax = (float) ((int) (maxValue / decimalFactor + 1.0))
				* decimalFactor;
		float scaleMin = (float) ((int) (minValue / decimalFactor))
				* decimalFactor;
		
		if(minValue < 0){
			scaleMin -= 1;
		}		

		return new float[] { scaleMin, scaleMax };
	}

	/**
	 * Calculates the number which is power of ten and just below preferredScale
	 * e.g. 1.3 -> 1, 23 -> 10, 543 -> 100
	 * 
	 * @param maxValue
	 * @return
	 */
	private float calculateDecimalFactor(float preferredScale) {
		// To scale the preset SCALES to the size of the prefferredScale
		float decimalFactor = 1;
		// radix is 10
		if (preferredScale >= 10) {
			while (preferredScale >= decimalFactor * 10) {
				decimalFactor *= 10;
			}
		} else {
			while (preferredScale < decimalFactor / 10) {
				decimalFactor /= 10;
			}
		}
		return decimalFactor;
	}

	public float[] getColorScaleValues() {
		return cScale;
	}

}
