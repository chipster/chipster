/*
 * Worker.java
 *
 * Created on 3. kesäkuuta 2006, 21:01
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Component;

import javax.swing.SwingUtilities;

/**
 * Class from the Viski project (http://www.cs.helsinki.fi/group/viski/). The newer Java
 * versions (from 1.5) might have the similar functionality already, but this isn't broken,
 * at least not too badly.
 *
 * @author esa
 */
public class Worker extends Thread {
    private Component paintable;
    private Projection projection;
    private boolean workRequested;
    private Object obj;
    
    
    /**
     * Creates a new instance of Worker
     * @param paintable 
     * @param projection 
     */
    public Worker(Component paintable, Projection projection) {
        this.paintable = paintable;
        this.projection = projection;
        this.workRequested = false;
        this.obj = new Object();
    }
    
    public void workRequest() {
        synchronized(obj) {
            workRequested = true;
            obj.notify();
        }
    }
    
    public void run() {
        while (true) {
            synchronized (obj) {
                while (workRequested == false) {
                    try {
                        obj.wait();
                    } catch (InterruptedException e) {}
                }
                
                SwingUtilities.invokeLater(new Runnable() {

                	@Override
                	public void run() {
                		projection.doProjection();
                		paintable.repaint();
                	}
                });                

                workRequested = false;
                obj.notify();
            }
        }
    }
}
