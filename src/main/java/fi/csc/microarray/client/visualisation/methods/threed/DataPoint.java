/*
 * DataPoint.java
 *
 * Created on 29. toukokuuta 2006, 23:58
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;

/**
 * Class for singel point of data. This is mainly from the Viski project 
 * (http://www.cs.helsinki.fi/group/viski/) but some new drawing modes are added.
 *
 * @author esa, Petri Klemelä
 */
public class DataPoint 
        extends Drawable {
    
    private double radius;
    private double effectiveRadius;
    private int index;
    
    private static final int SELECTION_DISTANCE = 4;
    
    /**
     * Creates a new instance of DataPoint
     * @param xData 
     * @param yData 
     * @param zData 
     * @param color 
     * @param radius 
     * @param text 
     */
    public DataPoint(double xData, double yData, double zData, 
            Color color, double radius, String text, int index) {
        this.radius = this.effectiveRadius = radius;
        this.visualisationCoords = new double[1][3];
        this.projectedCoords = new double[1][2];
        this.deviceCoords = new int[1][2];
        this.dataCoords = new double[1][3];
        this.distanceFromCamera = new double[1];
        this.index = index;        
        
        dataCoords[0][0] = visualisationCoords[0][0] = xData;
        dataCoords[0][1] = visualisationCoords[0][1] = yData;
        dataCoords[0][2] = visualisationCoords[0][2] = zData;
        this.color = color;
        this.text = text;
    }
    
    public int getIndex(){
    	return index;
    }
    
    /**
     * 
     * @param xData 
     * @param yData 
     * @param zData 
     * @param color 
     * @param radius 
     */
    public DataPoint(double xData, double yData, double zData, 
            Color color, double radius, CoordinateArea coordianteArea, String geneId) {
        this(xData, yData, zData, color, radius, null, -1);
    }

    /**
     * 
     * @param g 
     * @param width 
     * @param height 
     */
    public void draw(Graphics g, int width, int height, PaintMode paintMode) {
    	Graphics2D g2d = (Graphics2D)g; 
    	
        if (Math.abs(projectedCoords[0][0]) > 0.5 || 
                Math.abs(projectedCoords[0][1]) > 0.5 ||
                this.hidden == true) {
            return;
        }
        g2d.setColor(color);
        deviceCoords[0][0] = (int)((projectedCoords[0][0] + 0.5) * width);
        deviceCoords[0][1] = (int)((projectedCoords[0][1] + 0.5) * height);
        double screenRadius = 0.5;
        if (effectiveRadius * width > 0.5) {
            screenRadius = effectiveRadius * width;
        }       
        
        int x = deviceCoords[0][0]-(int)screenRadius;
        int y = deviceCoords[0][1]-(int)screenRadius;
        int w = (int)(screenRadius*2);
        int h = (int)(screenRadius*2);
     
        if(paintMode == CoordinateArea.PaintMode.RECT){
        	g2d.setPaint(color);        
        	g2d.fillRect(x, y, w, h);
        	
        } else {
        	paintBall(x, y, w, h, color, g2d);
        }
                              
        if (selected == true) {
        	g2d.setPaint(Color.gray);
        	g2d.drawOval(x-2, y-2, w+4, h+4);
        }
    }
    
    public static void paintBall(int x, int y, int w, int h, Color c, Graphics2D g2d) {    	  
    	
    	int radius = (int)(w * 1);
    	radius = radius > 0 ? radius : 1;
    	
    	g2d.setPaint(new RoundGradientPaint(
    		x + w / 4.0,
    		y + h / 4.0,
    		c,
    		new Point( radius, radius),
    		Color.BLACK));
    	
    	g2d.fillOval(x, y, w, h);
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
        
        double distanceScalar = planeDistance / distanceFromCamera[0];
        this.effectiveRadius = (radius / viewWindowWidth) * distanceScalar;
    }
    
    /**
     * 
     * @param x 
     * @param y 
     * @param points 
     * @param maxDist 
     * @return 
     */
    static final public LinkedList<DataPoint> getNearest(
            int x, 
            int y, 
            Drawable[] points, 
            double maxDist
            ) {
        LinkedList<DataPoint> list = new LinkedList<DataPoint>();
        DataPoint dp = null;
        double dist = Double.MAX_VALUE;
        
        for (Drawable d : points) {
            if (d instanceof DataPoint) {
                int xDiff = x-d.deviceCoords[0][0];
                int yDiff = y-d.deviceCoords[0][1];
                double tmp = Math.sqrt(xDiff*xDiff + yDiff*yDiff);
                if (tmp < dist && tmp <= maxDist) {
                    dist = tmp;
                    dp = (DataPoint)d;
                }
            }
        }
        if (dp != null)
            list.add(dp);
        return list;
    }
    
    /**
     * 
     * @param x 
     * @param y 
     * @param points 
     * @return 
     */
    static final public LinkedList<DataPoint> getNearest(
            int x, 
            int y, 
            Drawable[] points
            ) {
        return DataPoint.getNearest(x, y, points, SELECTION_DISTANCE);
    }
            
    /**
     * 
     * @param x1 
     * @param y1 
     * @param x2 
     * @param y2 
     * @param points 
     * @return 
     */
    static final public LinkedList<DataPoint> getGroup(
            int x1, 
            int y1, 
            int x2, 
            int y2, 
            Drawable[] points
            ) {
        LinkedList<DataPoint> list = new LinkedList<DataPoint>();
        //double dist = Double.MAX_VALUE;
        
        if (x1 > x2) {
            int tmp = x2;
            x2 = x1;
            x1 = tmp;
        }
        if (y1 > y2) {
            int tmp = y2;
            y2 = y1;
            y1 = tmp;
        }
        
        for (Drawable d : points) {
            if (d instanceof DataPoint) {
                if (d.deviceCoords[0][0] >= x1 && d.deviceCoords[0][0] <= x2 && 
                        d.deviceCoords[0][1] >= y1 && d.deviceCoords[0][1] <= y2) {
                    list.add((DataPoint)d);
                }
            }
        }
        return list;
    }
    
    /**
     * 
     * @param a 
     * @param b 
     * @return 
     */
    static final public double pointDistance(double[] a, double[] b) {
        double xD = a[0] - b[0];
        double yD = a[1] - b[1];
        double zD = a[2] - b[2];
        return Math.sqrt(xD*xD + yD*yD + zD*zD);
    }
    
    /*
     * In AnnotateListPanel this is used as gene identifier. If toString is changed, separate 
     * method for giving identifier should be done. 
     */
    @Override
    public String toString(){
    	return text;
    }
}
