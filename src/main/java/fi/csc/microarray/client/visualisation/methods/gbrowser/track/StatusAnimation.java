package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import javax.swing.JComponent;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * A Swing component which shows a spinning wheel when the background threads are working.
 * 
 * @author klemela
 */
public class StatusAnimation extends JComponent implements DataResultListener {
	
	private static final long VISIBLE_AFTER = 100; //ms
	private Map<DataThread, Long> queueLengths = new HashMap<DataThread, Long>();
	private double angle;
	private long previousTime;
	private long hideTime = 0;
	private QueueManager queueManager;
	private LinkedList<DataThread> dataThreads = new LinkedList<>();
	
	private static final int SIZE = 16;
	private static final int THICKNESS = 5;
	
	public StatusAnimation(QueueManager queueManager) {
		this.queueManager = queueManager;
	}

	@Override
	public void paintComponent(Graphics g) {
				
		long queueLength = getMaxQueueLength();
		
		//Show animation
		if (queueLength > 0) {
		
			if (System.currentTimeMillis() - hideTime > VISIBLE_AFTER) {
				for (double d = 0; d < 1; d += 0.1) {
					int radius = (SIZE - THICKNESS) / 2;
					double x = Math.sin(angle + d * Math.PI * 1.5) * radius;
					double y = -Math.cos(angle + d * Math.PI * 1.5) * radius;			

					angle += 0.005 * (System.currentTimeMillis() - previousTime);
					previousTime = System.currentTimeMillis();

					if (angle < 0) {
						angle += Math.PI * 2;
					}

					Color c = GBrowserConstants.COLOR_BLUE;
					c = new Color(c.getRed(), c.getGreen(), c.getBlue(), (int) (d * 128));

					Ellipse2D.Double circle = new Ellipse2D.Double(x + radius, y + radius,  THICKNESS, THICKNESS);
					g.setColor(c);
					((Graphics2D)g).fill(circle);
					
					//continue animation
					this.repaint();
				} 		
			}
		} else {			
			hideTime = System.currentTimeMillis();
		}				
	}
	
	@Override
	public Dimension getPreferredSize() {
		return new Dimension(SIZE, SIZE);
	}

	private long getMaxQueueLength() {
		long max = 0; 
		
		for (Long value : queueLengths.values()) {
			max = Math.max(max, value);
		}
		
		return max;
	}

	public void processDataResult(DataResult dataResult) {

		DataThread dataThread = dataResult.getStatus().getDataThread();
		
		if (dataResult.getStatus().getDataRequestCount() >= 0) {
						
			if (!queueLengths.containsKey(dataThread)) {
				queueLengths.put(dataThread, 0l);
			}

			Long value = dataResult.getStatus().getDataRequestCount(); 
			queueLengths.put(dataThread, value);
		}
	}

	public void addDataThread(DataThread dataThread) {
		this.dataThreads.add(dataThread);
	}

	public void initilizeListeners() {
		for (DataThread dataThread : dataThreads) {
			queueManager.addDataResultListener(dataThread, this);
		}	
	}

	public void clear() {
		dataThreads.clear();
	}
}
