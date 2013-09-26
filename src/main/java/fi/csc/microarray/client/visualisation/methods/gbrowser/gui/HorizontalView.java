package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.geom.Point2D;
import java.util.HashSet;
import java.util.Set;

import javax.swing.SwingUtilities;
import javax.swing.Timer;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;

/**
 * The basic view of genome browser. Shows the tracks horizontally.
 *
 */
public class HorizontalView extends GBrowserView implements KeyListener {

	private Timer keyTimer;
	private Set<Integer> keySet = new HashSet<Integer>();

	public HorizontalView(GBrowserPlot parent, boolean movable, boolean zoomable, boolean selectable) {
		super(parent, movable, zoomable, selectable);
		getComponent().addKeyListener(this);
		getComponent().addMouseListener(this);
		getComponent().addMouseMotionListener(this);
		getComponent().addMouseWheelListener(this);	
	}	

	@Override
	protected void handleDrag(Point2D start, Point2D end, boolean disableDrawing) {
		double bpMove = trackToBp((double) start.getX()).minus(trackToBp((double) end.getX()));

		RegionDouble newRegion = bpRegion.clone();
		
		newRegion.move(bpMove);
		setBpRegion(newRegion);

		if (!disableDrawing) {
			parentPlot.redraw();
		}
	}

	@Override
	public void keyPressed(KeyEvent e) {
			if ((zoomable || movable) && isNavigationKey(e.getKeyCode())) {
				keySet.add(e.getKeyCode());
				startKeyAnimation();				
			}
	}
	
	public boolean isNavigationKey(int keyCode) {
		return keyCode == KeyEvent.VK_UP ||
			keyCode == KeyEvent.VK_DOWN ||
			keyCode == KeyEvent.VK_LEFT ||
			keyCode == KeyEvent.VK_RIGHT;
	}
	
	public void startKeyAnimation() {
		
		if (keyTimer == null) {
			
			// Stop existing timers
			stopKeyAnimation();

			keyTimer = new Timer(1000 / FPS, new ActionListener() {

				//Frame indexes to keep speed constant regardless of performance
				private int i = 0;					

				private long startTime = System.currentTimeMillis();

				public void actionPerformed(ActionEvent arg0) {
					
					SwingUtilities.invokeLater(new Runnable() {
						public void run() {						

							// Always execute following loop at least once
							boolean skipFrame = true;
							
							// Last round has to draw
							boolean firstNotSkipped = true;	
							
							// Update locations without drawing until we are in the correct place for
							// moment of time and then do it once again with drawing
							while (skipFrame || firstNotSkipped){
																
								if (!skipFrame) {
									firstNotSkipped = false;
								}
																
								// Do until all keys are released or component has lost the focus
								if (!keySet.isEmpty() && getComponent().hasFocus()) {

									if (zoomable) {
										/* This value was obtained with trial-and-error method
										 * by trying to find value that would keep one side of the screen
										 * fixed when zooming and moving simultaneously. This could be 
										 * as well anything else as long as zooming speed feels nice.
										 */
										
										final double ZOOM_FACTOR = 0.78;
										
										if ( keySet.contains( KeyEvent.VK_UP )) {

											zoom(getWidth() / 2, -ZOOM_FACTOR, skipFrame);

										} 

										if ( keySet.contains( KeyEvent.VK_DOWN ))  {

											zoom(getWidth() / 2, ZOOM_FACTOR, skipFrame);
										}
									}

									if ( movable ) {
										
										final double SPEED_DIVIDER = 50.0;
										
										if ( keySet.contains( KeyEvent.VK_LEFT )) {
											bpRegion.move(-getBpRegion().getLength() / SPEED_DIVIDER);
											setBpRegion(bpRegion);
											
											if (!skipFrame) {
												parentPlot.redraw();											
											}
										} 

										if (keySet.contains(  KeyEvent.VK_RIGHT ))  {
											bpRegion.move(getBpRegion().getLength() / SPEED_DIVIDER);
											setBpRegion(bpRegion);

											if (!skipFrame) {
												parentPlot.redraw();											
											}
										}
									}

									i++;

								} else {

									stopKeyAnimation();
									break;
								}
								
								//Calculate if we are behind from our schedule and should just 
								//update the location without drawing on next round
								skipFrame = System.currentTimeMillis() > startTime + (1000.0 / FPS) * i;
							}
						}
					});
				}
			});

			keyTimer.setRepeats(true);
			keyTimer.setCoalesce(false);
			keyTimer.start();	
		}
	}

	private void stopKeyAnimation() {

		if (keyTimer != null) {
			keyTimer.stop();
			keyTimer = null;
			
			parentPlot.redraw();
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {
		
		if ((zoomable || movable) && isNavigationKey(e.getKeyCode())) {
			keySet.remove(e.getKeyCode());
		}
	}

	@Override
	public void keyTyped(KeyEvent e) {
		// ignore		
	}
}
