/*
 * Drawable.java
 *
 * Created on 30. toukokuuta 2006, 20:03
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Graphics;

import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;

/**
 * Super class for the items that are going to be drawn in Visualisation from the
 * Viski project (http://www.cs.helsinki.fi/group/viski/).
 *
 * @author esa
 */
abstract public class Drawable 
        implements Comparable<Drawable>{
    
    // Coordinates in camera space
    public double[][] visualisationCoords;
    // Coordinates on the projection plane
    public double[][] projectedCoords;
    // Coordinates on the screen
    public int[][] deviceCoords;
    // Coordinates in the original data space
    public double[][] dataCoords;
    
    public boolean selected = false;
    public boolean hidden = false;
    
    protected java.awt.Color color;
    protected String text;
    protected double[] distanceFromCamera;
    
    /**
     * 
     * @param g 
     * @param width 
     * @param height 
     * @param paintMode 
     */
    abstract public void draw(Graphics g, int width, int height, PaintMode paintMode);
    /**
     * 
     * @param camera 
     * @param planeDistance 
     * @param viewWindowWidth 
     * @param viewWindowHeight 
     */
    abstract public void setDistanceFromCamera(
            double[] camera, 
            double planeDistance,
            double viewWindowWidth,
            double viewWindowHeight);
    
    /**
     * 
     * @return 
     */
    public double[] getDistanceFromCamera() {
        return distanceFromCamera;
    }
    
    /**
	 * 
	 * @param other 
	 * @return 
	 */
	public int compareTo(Drawable other) {
		//TODO Find the reasen for the nullPointer    	
		if(other == null){
			//Avoid critical null pointer, which shouldn't really hapen
			System.out.println("null");
			return -1;
		}
	    double thisMax = Double.NEGATIVE_INFINITY;
	    double otherMin = Double.MAX_VALUE;
	    double[] d2 = other.getDistanceFromCamera();
	    
	    //Comapares maximum to minimum        
	    for (int i=0; i < distanceFromCamera.length; ++i) {
	        if (distanceFromCamera[i] > thisMax)
	            thisMax = distanceFromCamera[i];
	    }
	    for (int i=0; i < d2.length; ++i) {
	        if (d2[i] < otherMin)
	            otherMin = d2[i];
	    }
	    //Compares averages:
	    for (int i=0; i < distanceFromCamera.length; ++i) {            
	            thisMax += distanceFromCamera[i]/distanceFromCamera.length;
	    }
	    for (int i=0; i < d2.length; ++i) {
	            otherMin += d2[i]/d2.length;
	    }
	    
	    if (thisMax < otherMin)
	        return 1;
	    if (thisMax > otherMin)
	        return -1;
	    
	    return 0;
	}
}
