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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;

/**
 * From the Viski project (http://www.cs.helsinki.fi/group/viski/).
 *
 * @author esa
 */
public class Line 
        extends Drawable {
	
	private static final int COLOR_MARKER_WIDTH = 8;
	private Color textColor;
	private int thickness = 1;
	//To allow changing colors for the text    
	
    /**
     * Line with too ends and color
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
            double x2, double y2, double z2, Color color) {
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
    }
	
    /**
     * Line with too ends, text and color
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
            double x2, double y2, double z2, Color lineColor, String text, Color textColor) {
       
    	this(x1, y1, z1, x2, y2, z2, lineColor);        
        this.text = text;
        this.textColor = textColor;
    }

    public Line(double x1, double y1, double z1,
            double x2, double y2, double z2, Color lineColor, int thickness) {
    	this(x1, y1, z1, x2, y2, z2, lineColor);
    	this.thickness  = thickness;
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
        
        Graphics2D g2 = (Graphics2D)g;
        g2.setStroke(new BasicStroke(thickness));        
        g2.drawLine(deviceCoords[0][0], deviceCoords[0][1], 
                deviceCoords[1][0], deviceCoords[1][1]);
        if (this.text != null) {
        	
        	//Calculate text dimensions
        	int textWidth = 0;
        	int textHeight = 10;

     		textWidth += SwingUtilities.computeStringWidth(g.getFontMetrics(), text);
        	if (!Color.gray.equals(color)) {
        		textWidth += COLOR_MARKER_WIDTH;        		
        	}
        	
        	int textCenterX = deviceCoords[1][0];
        	int textCenterY = deviceCoords[1][1];        	        
        	
        	float dX = deviceCoords[1][0] - deviceCoords[0][0];
        	float dY = deviceCoords[1][1] - deviceCoords[0][1];
        	
        	float relativeX = textCenterX / (float)width; // form 0 to 1
        	float relativeY = textCenterY / (float)height; // from 0 to 1
        	        	
        	textCenterX += dX * Math.abs(relativeX - 0.5) * 1;
        	textCenterY += dY * relativeY * 0.2;        	

        	int x = (int) (textCenterX - textWidth / 2);
        	int y = textCenterY + textHeight;
        	
//        	int x = textCenterX - textWidth / 2;
//        	int y = textCenterY - textHeight / 2;
        	
        	//x = (int) (x - textWidth + textWidth * (relativeX - 0.25) * 2);
        	//y = (int) (y - textHeight + textHeight * (relativeY - 0.25) * 2);

        	int xShift = 0;        	        	

        	int textPartWidth = 0;
        	if (!Color.gray.equals(color)) {
        		g.setColor(textColor);
        		g.fillRect(x + xShift, y - 5, 6, 6);
        		textPartWidth += COLOR_MARKER_WIDTH;
        	}
        	g.setColor(Color.gray);
        	g.drawString(text, x + xShift + textPartWidth, y);
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
