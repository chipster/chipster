/*
 * Line.java
 *
 * Created on 30. toukokuuta 2006, 19:03
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;

/**
 * From the Viski project (http://www.cs.helsinki.fi/group/viski/).
 *
 * @author esa
 */
public class Line 
        extends Drawable {
	
	//To allow changing colors for the text
	private Color[] textColors;
	private String[] texts;
    
    /**
     * Creates a new instance of Line
     * @param x1 
     * @param y1 
     * @param z1 
     * @param x2 
     * @param y2 
     * @param z2 
     * @param color 
     * @param text 
     */
	
	
    /**
     * Line with too ends, text and same color for all these
     * 
     * @param x1
     * @param y1
     * @param z1
     * @param x2
     * @param y2
     * @param z2
     * @param color
     * @param text
     */
    public Line(double x1, double y1, double z1,
            double x2, double y2, double z2, Color color, String text) {
        this.visualisationCoords = new double[2][3];
        this.projectedCoords = new double[2][2];
        this.deviceCoords = new int[2][2];
        this.dataCoords = new double[2][3];
        this.distanceFromCamera = new double[2];
        
        dataCoords[0][0] = visualisationCoords[0][0] = x1;
        dataCoords[0][1] = visualisationCoords[0][1] = y1;
        dataCoords[0][2] = visualisationCoords[0][2] = z1;
        dataCoords[1][0] = visualisationCoords[1][0] = x2;
        dataCoords[1][1] = visualisationCoords[1][1] = y2;
        dataCoords[1][2] = visualisationCoords[1][2] = z2;
        this.color = color;
        
        this.texts = new String[]{ text };
        this.textColors = new Color[]{ color };
    }
    
	
    /**    
     * Line with too ends, line color, and text with multiple colors
     * 
     * @param x1
     * @param y1
     * @param z1
     * @param x2
     * @param y2
     * @param z2
     * @param color
     * @param texts Array of text which is shown concatenated
     * @param textColors Array of colors, indexes correspond to the indexes in array texts.
     */
    public Line(double x1, double y1, double z1,
            double x2, double y2, double z2, 
            Color color, String[] texts, Color[] textColors){
    	
    	this(x1, y1, z1, x2, y2, z2, color, null);
    	this.texts = texts;
    	this.textColors = textColors;
    }
    
    /**
     * 
     * @param g 
     * @param width 
     * @param height 
     */
    public void draw(Graphics g, int width, int height, PaintMode notUsed) {
    	
    	//Disable drawing lines very near camera, because it causes strange effects.
    	//Value of 1.5 was found with visual experiment, so feel free to modify it.
        if (Math.abs(projectedCoords[0][0]) > 1.5 || 
                Math.abs(projectedCoords[0][1]) > 1.5 ||
                this.hidden) {
            return;
        }
        
        g.setColor(color);
        deviceCoords[0][0] = (int)((projectedCoords[0][0] + 0.5) * width);
        deviceCoords[0][1] = (int)((projectedCoords[0][1] + 0.5) * height);
        deviceCoords[1][0] = (int)((projectedCoords[1][0] + 0.5) * width);
        deviceCoords[1][1] = (int)((projectedCoords[1][1] + 0.5) * height);
        
        g.drawLine(deviceCoords[0][0], deviceCoords[0][1], 
                deviceCoords[1][0], deviceCoords[1][1]);
        if (this.texts != null) {
        	
        	//Calculate text dimension
        	int textWidth = 0;
        	int textHeight = 10;
        	for(String text : texts){
        		textWidth += SwingUtilities.computeStringWidth(g.getFontMetrics(), text); 
        	}
        	
        	//Place text from the line end away from the center of the view area
        	textWidth *=  (1.5 -  2*deviceCoords[1][0] / (float)width);
        	
        	if(deviceCoords[1][1] < height / 2.0){
        		textHeight = 0;
        	}        
        	
        	int x = deviceCoords[1][0] - textWidth;
        	int y = deviceCoords[1][1] + textHeight -1;
        	
        	int xShift = 0;
        	for(int i = 0; i < texts.length && i < textColors.length ; i++){
        		g.setFont(g.getFont().deriveFont(9.0f));        	        		
        		g.setColor(textColors[i]);
        		g.drawString(texts[i], x + xShift, y);        		
        		xShift += SwingUtilities.computeStringWidth(g.getFontMetrics(), texts[i]);
        	}
        }
    }

    /**
     * 
     * @param camera 
     * @param planeDistance 
     * @param viewWindowWidth 
     * @param viewWindowHeight 
     */
    public void setDistanceFromCamera(
            double[] camera,
            double planeDistance,
            double viewWindowWidth,
            double viewWindowHeight) {
        this.distanceFromCamera[0] = 
                DataPoint.pointDistance(camera, visualisationCoords[0]);
        this.distanceFromCamera[1] = 
                DataPoint.pointDistance(camera, visualisationCoords[1]);
    }

}
