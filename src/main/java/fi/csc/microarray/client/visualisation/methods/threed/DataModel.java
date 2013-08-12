package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.util.ScaleUtil;

/**
 * Class gets the data values from the iterables and scales them to the scales that
 * are also decided here.
 * 
 * @author Petri Klemelä
 */
public class DataModel {
	
	private static final Logger logger = Logger
	.getLogger(DataModel.class);
	
	private Color lineColor = Color.DARK_GRAY;
	
	// Color schemes from http://colorbrewer2.org/  Apache License, Version 2.0	
	//11-class Paired qualitative
	public static Color[] qualitatativeColorScheme = new Color[] {
		
		//custom order: intense colors first
		new Color(31, 120, 180), 
		new Color(51, 160, 44),
		new Color(227, 26, 28),
		new Color(255, 127, 0),
		new Color(106, 61, 154),
		new Color(166, 206, 227), 
		new Color(178, 223, 138), 
		new Color(251, 154, 153), 
		new Color(253, 191, 111), 
		new Color(202, 178, 214), 
		new Color(255, 255, 153),
		
		Color.lightGray,
		Color.darkGray,
	};
		
//		original order
//		new Color(166, 206, 227), 
//		new Color(31, 120, 180), 
//		new Color(178, 223, 138), 
//		new Color(51, 160, 44),
//		new Color(251, 154, 153), 
//		new Color(227, 26, 28),
//		new Color(253, 191, 111), 
//		new Color(255, 127, 0),
//		new Color(202, 178, 214), 
//		new Color(106, 61, 154),
//		new Color(255, 255, 153) 
//	};

	// 5-class Yellow-Green-Blue sequential
//	Color[] sequentialColorScheme = new Color[] {
//			new Color(255, 255, 204), 
//			new Color(161, 218, 180), 
//			new Color(65, 182, 196), 
//			new Color(44, 127, 184), 
//			new Color(37, 52, 148),
//	};	
	
	//5-class Red-Purple sequential
//	Color[] sequentialColorScheme = new Color[] {
//		new Color(254, 235, 226), 
//		new Color(251, 180, 185),
//		new Color(247, 104, 161),
//		new Color(197, 27, 138),
//		new Color(122, 1, 119)
//	};
	
	//9-class Yellow-Orange-Red
	public static Color[] sequentialColorScheme = new Color[] {
			new Color(255, 255, 204), 	
			new Color(255, 237, 160), 
			new Color(254, 217, 118),
			new Color(254, 178, 76),
			new Color(253, 141, 60),
			new Color(252, 78, 42),
			new Color(227, 26, 28),
			new Color(189, 0, 38),
			new Color(128, 0, 38)
	};

	/**
	 * @param value
	 *            between 0 and 1.0
	 * @return corresponding color
	 */
	public Color getColorFor(float value) {
		Color c = Color.gray;
		for (int i = 0; i < cScale.length - 1; i++) {
			if (value >= cScale[i]) {
				c = getColorScheme()[i];
			}
		}
		return c;
	}

	private Drawable[] points;

	//How many lines are shown in the grid for each axis
	private static final int LINE_COUNT = ScaleUtil.DEFAULT_STEP_COUNT;

	private float[] cScale;

	private Float[] cValues;
	private Color[] colorScheme;
	
	public DataModel(){			
		this(sequentialColorScheme);
	}
		
	
	public DataModel(Color[] colorScheme) {
		this.colorScheme = colorScheme;
	}

	public void setData(Iterable<String> identifiers, Iterable<Float> xValues, Iterable<Float> yValues,
			Iterable<Float> zValues, Iterable<Float> cValuesIterable) {			
		
		List<Drawable> points = new ArrayList<Drawable>();

		Iterator<String> identIter = identifiers.iterator();
		Iterator<Float> xIter = xValues.iterator();
		Iterator<Float> yIter = yValues.iterator();
		Iterator<Float> zIter = zValues.iterator();
		Iterator<Float> cIter = cValuesIterable.iterator();
		
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

		float[] xScale = ScaleUtil.generateScaleValues(minMaxX[0], minMaxX[1]);
		float[] yScale = ScaleUtil.generateScaleValues(minMaxY[0], minMaxY[1]);
		float[] zScale = ScaleUtil.generateScaleValues(minMaxZ[0], minMaxZ[1]);
		cScale = this.generateColorScaleValues(minMaxC[0], minMaxC[1]);

		points.addAll(this.generateScaleLines(xScale, yScale, zScale));

		// reset iterators
		xIter = xValues.iterator();
		yIter = yValues.iterator();
		zIter = zValues.iterator();
		cIter = cValuesIterable.iterator();

		logger.debug((System.currentTimeMillis() - time) / 1000);
		time = System.currentTimeMillis();
		
		//TODO very slow, 25 secs with 22000 rows
		xIter = xValues.iterator();
		yIter = yValues.iterator();
		zIter = zValues.iterator();
		cIter = cValuesIterable.iterator();
		
		List<Float> cList = new LinkedList<Float>();
		
		for (int i = 0; identIter.hasNext() && xIter.hasNext() && yIter.hasNext() && 
				zIter.hasNext() && cIter.hasNext(); i++) {
			String identifier = identIter.next();
			float x = xIter.next();
			float y = yIter.next();
			float z = zIter.next();
			float c = cIter.next();
			float[] scaled = this.convertToScaled(xScale, yScale, zScale,
					cScale, new float[] { x, y, z, c });
			
			cList.add(c);
			
			Color color = getColorFor(scaled[3]);

			points.add(new DataPoint(scaled[0], scaled[1], scaled[2], color,
					0.015, identifier, i));
		}
		
		this.cValues = cList.toArray(new Float[0]);
				
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
		float max = Float.NEGATIVE_INFINITY;
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
				original[3]};
	}
	
	protected float convertToScaled(float scale[], float value){
		return (value - scale[0]) / (scale[LINE_COUNT - 1] - scale[0]);
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
	private String coordToString(float x, float y, float z){
		return 
				"(" +  
				"x " + ScaleUtil.format(x) +  
				", " + 
				"y " + ScaleUtil.format(y) +  
				", " + 
				"z " + ScaleUtil.format(z) +  
				")"; 
	}

	private List<Drawable> generateScaleLines(float[] xValues, float[] yValues,
			float[] zValues) {	
				
		Color xColor = Color.gray;
		Color yColor = Color.gray;
		Color zColor = Color.gray;
		Color gray = Color.gray;
		
		//To understand following  method calls it's beneficial to draw
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

		//Coordinate axes
		List<Drawable> drawables = new ArrayList<Drawable>();
		drawables.add(new Line(0, 0, 0, 1, 0, 0, lineColor, 3));				
		drawables.add(new Line(0, 0, 0, 0, 1, 0, lineColor, 3));
		drawables.add(new Line(0, 0, 0, 0, 0, 1, lineColor, 3));
		
		//arrows of the coordinate axes
		drawables.add(new Line(1,0,0, 0.99, 0.005, -0.01, lineColor, 3));
		drawables.add(new Line(1,0,0, 0.99, -0.01, 0.005, lineColor, 3));
		
		drawables.add(new Line(0,1,0, 0.005, 0.99, -0.01, lineColor, 3));
		drawables.add(new Line(0,1,0, -0.01, 0.99, 0.005, lineColor, 3));
		
		drawables.add(new Line(0,0,1, -0.01, 0.005, 0.99, lineColor, 3));
		drawables.add(new Line(0,0,1, 0.005, -0.01, 0.99, lineColor, 3));
		
		//labels of the coordinate axes
		drawables.add(new Text(1.03f, 0, 0, "X", gray, true));				
		drawables.add(new Text(0, 1.03f, 0, "Y", gray, true));
		drawables.add(new Text(0, 0, 1.03f, "Z", gray, true));
		

		//Cube corner labels
		final float TEXT_0 = 0f;
		final float TEXT_1 = 1.1f;
		
		int L = LINE_COUNT - 1; // last index
		
		drawables.add(new Text(TEXT_0, TEXT_0, TEXT_0, coordToString(xValues[0], yValues[0], zValues[0]), gray));
		drawables.add(new Text(TEXT_1, TEXT_0, TEXT_0, coordToString(xValues[L], yValues[0], zValues[0]), gray));
		drawables.add(new Text(TEXT_1, TEXT_0, TEXT_1, coordToString(xValues[L], yValues[0], zValues[L]), gray));
		drawables.add(new Text(TEXT_0, TEXT_0, TEXT_1, coordToString(xValues[0], yValues[0], zValues[L]), gray));
		
		drawables.add(new Text(TEXT_0, TEXT_1, TEXT_0, coordToString(xValues[0], yValues[L], zValues[L]), gray));
		drawables.add(new Text(TEXT_1, TEXT_1, TEXT_0, coordToString(xValues[L], yValues[L], zValues[0]), gray));
		//lines.add(new Text(TEXT_1, TEXT_1, TEXT_1, coordToStringTable(xValues[L], yValues[L], zValues[L]), gray));
		drawables.add(new Text(TEXT_0, TEXT_1, TEXT_1, coordToString(xValues[0], yValues[L], zValues[L]), gray));

		//Cube side labels
		String xText;
		String yText;
		String zText;

		for (int i = 1; i < LINE_COUNT; i++) {
			
			float pos = (float)i / (LINE_COUNT - 1); 

			if (i != LINE_COUNT - 1) {
				//Normal lines show only one coordinate
				xText = "x " + ScaleUtil.format(xValues[i]);
				yText = "y " + ScaleUtil.format(yValues[i]);
				zText = "z " + ScaleUtil.format(zValues[i]);
				
			} else {
				//Mark texts disabled ( the last lines )
				xText = null;
				yText = null;
				zText = null;								
			}

			// X-scale 
			drawables.add(new Line(pos, 0, 0, pos, 1, 0, lineColor));
			drawables.add(new Line(pos, 0, 0, pos, 0, 1, lineColor));
			
			drawables.add(new Text(pos, TEXT_1, TEXT_0, xText, xColor));
			drawables.add(new Text(pos, TEXT_0, TEXT_1, xText, xColor));


			// Y-scale 
			drawables.add(new Line(0, pos, 0, 0, pos, 1, lineColor));
			drawables.add(new Line(0, pos, 0, 1, pos, 0, lineColor));
			
			drawables.add(new Text(TEXT_0, pos, TEXT_1, yText, yColor));
			drawables.add(new Text(TEXT_1, pos, TEXT_0, yText, yColor));

			// Z-Scale 
			drawables.add(new Line(0, 0, pos, 1, 0, pos, lineColor));
			drawables.add(new Line(0, 0, pos, 0, 1, pos, lineColor));
			
			drawables.add(new Text(TEXT_1, TEXT_0, pos, zText, zColor));
			drawables.add(new Text(TEXT_0, TEXT_1, pos, zText, zColor));

		}

		return drawables;
	}
	
	private float[] generateColorScaleValues(float min, float max) {
		
		int colorCount = getColorScheme().length + 1;
		
		float[] values = new float[colorCount];
		for (int i = 0; i < colorCount; i++) {
			values[i] = min + (max - min) / (colorCount - 1) * i;
		}

		return values;
	}


	public float[] getColorScaleValues() {
		return cScale;
	}
	
	public Float[] getColorValues() {
		return cValues;
	}
	
	public Color[] getColorScheme() {
		return colorScheme;
	}
}
