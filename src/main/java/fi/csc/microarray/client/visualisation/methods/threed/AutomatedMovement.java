/*
 * AutomatedMovement.java
 *
 * Created on 3. kesäkuuta 2006, 20:53
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.LinkedList;
import java.util.Random;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;

/**
 * Thiss class takes care of the rotation animations. Base functionality is from the Viski library 
 * (http://www.cs.helsinki.fi/group/viski/) but lot of refactoring and bugfixes are done. Also
 * support for kinetic movement is added.
 *
 * @author esa, Petri Klemelä
 */
public class AutomatedMovement extends Thread {

    private Projection projection;
    private LinkedList<Task> taskQueue;
    //private Object obj;
    
    
    
    RotationTask rotationTask;
    Task task;
	private boolean kill;
	private CoordinateArea coordinateArea;
    
    /**
     * Creates a new instance of AutomatedMovement
     * @param projection 
     * @param coordinateArea 
     */
    public AutomatedMovement(Projection projection, CoordinateArea coordinateArea) {    	    	
        this.projection = projection;
        this.taskQueue = new LinkedList<Task>();
        this.coordinateArea = coordinateArea;
        //this.obj = new Object();
        
        Session.getSession().getApplication().addClientEventListener(new PropertyChangeListener(){
			public void propertyChange(PropertyChangeEvent evt) {
				if(evt instanceof VisualisationMethodChangedEvent){
					AutomatedMovement.this.kill();
				}
			}        	
        });
    }
    
    /**
     * 
     * @param point 
     * @param time 
     * @param fps 
     */
    public void addTranslationTask (double[] point, long time, double fps) {
        synchronized(taskQueue) {
            taskQueue.addLast(new TranslationTask(point, time, fps));
            taskQueue.notify();
        }
    }
    
    public void startAutomatedRotation(){
    	synchronized(taskQueue) {
        	
        	//Stop the previous tasks
    		clearTasks();
        	
            taskQueue.addLast(new AutomatedRotation());
            taskQueue.notify();
        }
    }        
    
    /**
     * 
     * @param xRotation 
     * @param yRotation 
     * @param zRotation 
     * @param time 
     * @param fps 
     */
    public void addRotationTask (double xRotation, double yRotation, double zRotation, 
            long time, double fps) {
        synchronized(taskQueue) {
        	
        	//Stop the previus tasks
        	clearTasks();
        	
            taskQueue.addLast(new RotationTask(xRotation, yRotation, zRotation, time, fps));
            taskQueue.notify();
        }
    }
    
    public RotationTask startKineticMove(double fps, double retardation){
    	
    	synchronized(taskQueue) {

    		//Stop the previus tasks
    		clearTasks();
    		rotationTask = new RotationTask(0, 0, 0, fps, retardation);
    		taskQueue.addLast(rotationTask);
    		taskQueue.notify();
    	}
    	return rotationTask;
    }
    
    public void clearTasks(){
    	while(taskQueue.size() > 0){
			taskQueue.remove();
		}
    	
    	if(task != null){
    		task.stop();
    	}
    }
    
    public RotationTask restartKineticMove(){
    	
    	synchronized(taskQueue) {

    		//Stop the previus tasks
    		while(taskQueue.size() > 0){
    			taskQueue.remove();
    		}

    		taskQueue.addLast(rotationTask);
    		taskQueue.notify();
    	}
    	return rotationTask;
    }
    
    public void run() {
        
        while (!kill) {
            synchronized(taskQueue) {
                while (taskQueue.isEmpty()) {
                    try {
                    	taskQueue.wait();
                    } catch (InterruptedException e) {}
                }
                task = taskQueue.poll();
                taskQueue.notify();
            }
            
            task.doTask();
            
            task = null;
        }
    }
    
    private abstract class Task {
        public abstract void doTask();        
        public abstract void stop();
    }
    
    private class TranslationTask extends Task {
        private double[] point;
        private long time;
        private double fps;
        
        public TranslationTask(double[] point, long time, double fps) {
            this.point = point;
            this.time = time;
            this.fps = fps;
        }
        
        public void stop(){
        	//TODO implement if transtlation is taken in use
        }
        
        public void doTask() {
            double ticks = (time/1000.0)*fps;
            double timeIncrement = 1000.0/fps;
            double[] placeIncrement = new double[3];
            placeIncrement[0] = point[0]/ticks;
            placeIncrement[1] = point[1]/ticks;
            placeIncrement[2] = point[2]/ticks;
            
            for (int i=0; i < ticks; ++i) {
                projection.setNewOrigin(placeIncrement);
                repaint();
                try {
                    sleep((long)timeIncrement);
                }
                catch (InterruptedException e) {}
            }
           
        }
        
    }
    
    public class RotationTask extends Task {
                
		double ticks;
        double timeIncrement;
        double xAngleInc;
        double yAngleInc;
        double zAngleInc;   
        double retardation = 1;
        
        
        
        //Task with target
        public RotationTask(double xRotation, double yRotation, double zRotation, 
                long time, double fps) {
            
            ticks = (time/1000.0)*fps;
            timeIncrement = 1000.0/fps;
            xAngleInc = (xRotation-projection.getXAxisRotation())/ticks;
            yAngleInc = (yRotation-projection.getYAxisRotation())/ticks;
            zAngleInc = (zRotation-projection.getZAxisRotation())/ticks;   
        }
        
        //Task with initial speed
        public RotationTask(double xAngleInc, double yAngleInc, double zAngleInc, double fps, double retardation) {
            
            ticks = Long.MAX_VALUE;            
            timeIncrement = 1000.0/fps;
            this.retardation = retardation;
            
            this.xAngleInc = xAngleInc;
            this.yAngleInc = yAngleInc;
            this.zAngleInc = zAngleInc;             
        }        
        
        public void setAngleIncs(double xAngleInc, double yAngleInc, double zAngleInc) {

        	ticks = Long.MAX_VALUE;
        	
            this.xAngleInc = xAngleInc;
            this.yAngleInc = yAngleInc;
            this.zAngleInc = zAngleInc;             
        }
        
        public void stop(){
        	ticks = -1;
        }
        
        public void doTask() {                
            
            for (int i=0; i < ticks; ++i) {
                                  	
                projection.setXAxisRotation(projection.getXAxisRotation()+xAngleInc);
                projection.setYAxisRotation(projection.getYAxisRotation()+yAngleInc);
                projection.setZAxisRotation(projection.getZAxisRotation()+zAngleInc);
                repaint();
                
                this.xAngleInc *= retardation;
                this.yAngleInc *= retardation;
                this.zAngleInc *= retardation;
                
                try {
                    sleep((long)timeIncrement);
                }
                catch (InterruptedException e) {}
            }
        }        
    }
    
	private class AutomatedRotation extends Task {
		
		boolean autoRotationKill = false;
	    boolean autoRotationIsRunning = false;
		
		double x = projection.getXAxisRotation();
		double y = projection.getYAxisRotation();
		double z = projection.getZAxisRotation();
		
		double angle = Math.PI*3/2;
		final double SPEED = 0.01;
		Random rand = new Random();		
		
		public void doTask(){
			
			autoRotationIsRunning = true;
			while(!autoRotationKill){
				
				projection.setXAxisRotation(x);
				projection.setYAxisRotation(y);
				projection.setZAxisRotation(z);

				x += Math.cos(angle) * SPEED;
				y += Math.sin(angle) * SPEED;
				angle += (rand.nextDouble()-0.5)*0.1;

				repaint();

				try {
					Thread.sleep(50);
				} catch (InterruptedException e) {
				}
			}			
			autoRotationIsRunning = false;
			autoRotationKill = false;
		}
		
		public void stop(){
			if(autoRotationIsRunning){
				this.autoRotationKill = true;
			}	
		}
	}
	
	private void repaint() {
		SwingUtilities.invokeLater(new Runnable(){
			public void run(){
				coordinateArea.repaint();
			}
		});
	}

	public void kill() {
		this.clearTasks();
		kill = true;
	}    
	
}
