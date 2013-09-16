/*
 * Projection.java
 *
 * Created on 28. toukokuuta 2006, 19:59
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.util.concurrent.PriorityBlockingQueue;

/**
 * Takes care of the current projection parameters, as rotation and camera
 *
 * @author esa, Petri Klemelä
 */
public class Projection {
    DataModel dataModel;
    /** Creates a new instance of Projection */
    public Projection(DataModel dataModel) {
    	this.dataModel = dataModel;
    }
    
    private double xAxisRotation = Math.PI*1;
    private double yAxisRotation = Math.PI*0;
    private double zAxisRotation = Math.PI*0;
    
    //Angle of view is correlative of these
    private double distanceOfProjectionPlaneFromCamera = 15;
    private double viewWindowWidth = 5;
    private double viewWindowHeight = 5;
    
    private double[] camera = {0,0,-5};
    private double[] newOrigin = {0.5,0.5,0.5};
    
    private PriorityBlockingQueue<Drawable> resultPoints = null;   
    
    /**
     * 
     * @return 
     */
    public double getXAxisRotation() {
        return xAxisRotation;
    }

    /**
     * 
     * @param xAxisRotation 
     */
    public void setXAxisRotation(double xAxisRotation) {
        this.xAxisRotation = xAxisRotation;
    }

    /**
     * 
     * @return 
     */
    public double getYAxisRotation() {
        return yAxisRotation;
    }

    /**
     * 
     * @param yAxisRotation 
     */
    public void setYAxisRotation(double yAxisRotation) {
        this.yAxisRotation = yAxisRotation;
    }

    /**
     * 
     * @return 
     */
    public double getZAxisRotation() {
        return zAxisRotation;
    }

    /**
     * 
     * @param zAxisRotation 
     */
    public void setZAxisRotation(double zAxisRotation) {
        this.zAxisRotation = zAxisRotation;
    }

    /**
     * 
     * @return 
     */
    public double getDistanceOfProjectionPlaneFromOrigin() {
        return distanceOfProjectionPlaneFromCamera;
    }

    /**
     * 
     * @param distanceOfProjectionPlaneFromOrigin 
     */
    public void setDistanceOfProjectionPlaneFromOrigin(
            double distanceOfProjectionPlaneFromOrigin
            ) {
        this.distanceOfProjectionPlaneFromCamera = 
                distanceOfProjectionPlaneFromOrigin;
    }

    /**
     * 
     * @return 
     */
    public double getViewWindowWidth() {
        return viewWindowWidth;
    }

    /**
     * 
     * @param viewWindowWidth 
     */
    public void setViewWindowWidth(double viewWindowWidth) {
        this.viewWindowWidth = viewWindowWidth;
    }

    /**
     * 
     * @return 
     */
    public double getViewWindowHeight() {
        return viewWindowHeight;
    }

    /**
     * 
     * @param viewWindowHeight 
     */
    public void setViewWindowHeight(double viewWindowHeight) {
        this.viewWindowHeight = viewWindowHeight;
    }

    /**
     * 
     * @return 
     */
    public double[] getPointOfView() {
        return camera;
    }

    /**
     * 
     * @param pointOfView 
     */
    public void setPointOfView(double[] pointOfView) {
        this.camera = pointOfView;
    }
    
    /**
     * 
     * @return 
     */
    public double[] getNewOrigin() {
        return newOrigin;
    }

    /**
     * 
     * @param newOrigin 
     */
    public void setNewOrigin(double[] newOrigin) {
        this.newOrigin[0] = newOrigin[0];
        this.newOrigin[1] = newOrigin[1];
        this.newOrigin[2] = newOrigin[2];
    }
    
    /**
     * 
     * @return 
     */
    public PriorityBlockingQueue<Drawable> getResultPoints() {
        return resultPoints;
    }

    private Matrix xAxisRotationMatrix(double degree) {
        double cosDegree = Math.cos(degree);
        double sinDegree = Math.sin(degree);
        Matrix m = new Matrix(4,4,0);
        
        m.set(1,1,cosDegree);
        m.set(2,2,cosDegree);
        m.set(1,2,-sinDegree);
        m.set(2,1,sinDegree);
        m.set(0,0,1);
        m.set(3,3,1);
        
        return m;
    }
    
    private Matrix yAxisRotationMatrix(double degree) {
        double cosDegree = Math.cos(degree);
        double sinDegree = Math.sin(degree);
        Matrix m = new Matrix(4,4,0);
        
        m.set(0,0,cosDegree);
        m.set(2,2,cosDegree);
        m.set(2,0,-sinDegree);
        m.set(0,2,sinDegree);
        m.set(1,1,1);
        m.set(3,3,1);
        
        return m;
    }
    
    private Matrix zAxisRotationMatrix(double degree) {
        double cosDegree = Math.cos(degree);
        double sinDegree = Math.sin(degree);
        Matrix m = new Matrix(4,4,0);
        
        m.set(0,0,cosDegree);
        m.set(1,1,cosDegree);
        m.set(0,1,-sinDegree);
        m.set(1,0,sinDegree);
        m.set(2,2,1);
        m.set(3,3,1);
        
        return m;
    }
    
    private Matrix translatePointToOriginMatrix(double[] point) {
        Matrix m = new Matrix(4,4,0);
        
        m.set(0,0,1);
        m.set(1,1,1);
        m.set(2,2,1);
        m.set(3,3,1);
        m.set(0,3,-point[0]);
        m.set(1,3,-point[1]);
        m.set(2,3,-point[2]);
        
        return m;
    }
    
    private Matrix perspectiveProjectionMatrix(double planeDistance) {
        Matrix m = new Matrix(4,4,0);
        
        m.set(0,0,1);
        m.set(1,1,1);
        m.set(2,2,1);
        m.set(3,2,1/planeDistance);
        
        return m;
    }
    
    private void homogeneousPointToProjectionPlanePoint(Matrix m, double[] p) {
        double w = m.get(3,0);
        p[0] = m.get(0,0)/w;
        p[1] = m.get(1,0)/w;
        //point3D[2] = m.get(2,0)/w;
    }
    
    /**
     * 
     * @return 
     */
    public PriorityBlockingQueue<Drawable> doProjection() {
        Matrix projM = perspectiveProjectionMatrix(distanceOfProjectionPlaneFromCamera);
                
        projM = projM.times(translatePointToOriginMatrix(camera));
        
        //Order of these decides which rotation is done around camera axis and
        //which around initial axis. 
        Matrix rotM = zAxisRotationMatrix(zAxisRotation);        
        rotM = rotM.times(xAxisRotationMatrix(xAxisRotation));
        rotM = rotM.times(yAxisRotationMatrix(yAxisRotation));
                
        rotM = rotM.times(translatePointToOriginMatrix(newOrigin));
        
        Drawable[] points = dataModel.getDataArray();
        
        resultPoints = new PriorityBlockingQueue<Drawable>(points.length);
 
        for (int i=0; i < points.length; i++) {
            boolean visible = rotateAndProject(points[i], rotM, projM, camera,
                    distanceOfProjectionPlaneFromCamera,
                    viewWindowWidth, viewWindowHeight);
            
            if (visible == true && points[i] != null) {
                resultPoints.add(points[i]);
            }
        }
        
        return resultPoints;
    }
    
    private boolean rotateAndProject(Drawable d, Matrix rot, Matrix proj,
            double[] cam, double planeDist, double viewW, double viewH) {
        Matrix point = new Matrix(4,1,1);
        
        boolean visible = false;
        
        for (int i=0; i < d.visualisationCoords.length; ++i) {
        	
        	//To make rotations not cumulative
        	d.visualisationCoords[i][0] = d.dataCoords[i][0];
        	d.visualisationCoords[i][1] = d.dataCoords[i][1];
        	d.visualisationCoords[i][2] = d.dataCoords[i][2];
        	
            point.set(0,0,d.visualisationCoords[i][0]);
            point.set(1,0,d.visualisationCoords[i][1]);
            point.set(2,0,d.visualisationCoords[i][2]);
            
            point = rot.times(point);
            
            //To get size of the ball right
            d.visualisationCoords[i][0] = point.get(0,0);
            d.visualisationCoords[i][1] = point.get(1,0);
            d.visualisationCoords[i][2] = point.get(2,0);
            
            
            point = proj.times(point);
            
            homogeneousPointToProjectionPlanePoint(point, d.projectedCoords[i]);
            
            d.projectedCoords[i][0] /= viewW;
            d.projectedCoords[i][1] /= viewH;
            
             if (d.visualisationCoords[i][2] > camera[2] /*+ planeDist ???*/)
                 visible = true;
            
            point.set(3, 0, 1);
        }
        
        d.setDistanceFromCamera(cam, planeDist,viewW, viewH);
        return visible;
    }

}
