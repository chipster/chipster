package fi.csc.microarray.client.dataview;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Timer;

import org.apache.log4j.Logger;


/**
 * Class to make automated movements of the graph smoother
 * 
 * @author klemela
 */
public class ScrollAnimator implements ActionListener{
	private Point from;
	private Point viewTarget = new Point();
	AnimatorScrollable scrollable;
	
	//to avoid dividing by zero
	private int frame = 1;
	private Timer timer;
	private static boolean isBusy = false;
	private static ScrollAnimator latestAnimator;

	public static final int FRAME_COUNT = 10;
	public static final int SCROLL_MARGIN = 10;
	
	private static final Logger logger = Logger.getLogger(ScrollAnimator.class);

	/**
	 * @param graphPanel to get graph and scroller objects
	 * @param to - Bounds of the component we want to show
	 */
	public ScrollAnimator(GraphPanel graphPanel, Rectangle to){
		this((Point)graphPanel.getGraph().fromScreen(graphPanel.getScroller().getViewport().getViewPosition()),//From
				to,
				graphPanel.getGraph().fromScreenCoordinates(graphPanel.getScroller().getViewport().getSize()), //viewSize
				graphPanel //scrollable
				);
	}
		
	public ScrollAnimator(Point from, Rectangle to, Dimension viewSize, AnimatorScrollable scrollable){				
		logger.debug("Animator started");
		if(latestAnimator != null){
			latestAnimator.kill();
		}
		latestAnimator = this;
		
		this.scrollable = scrollable; 
		this.from = from;		
		
		//If the target is visible, don't do anything
		Rectangle area = new Rectangle(from,viewSize);
		//area.grow(-SCROLL_MARGIN, -SCROLL_MARGIN);
		if(area.contains(to)){
			logger.debug("Showing allready");
			return;
		}
		
		//Depending on direction of scroll, move minimum amount required to show 
		//the target area
		if(from.getX() < to.getX()){
			viewTarget.x = (int)(to.getX()+to.getWidth()-viewSize.getWidth()+
					SCROLL_MARGIN);
		}else{
			viewTarget.x = (int)(to.getX()-SCROLL_MARGIN);
		}
		
		if(from.getY() < to.getY()){
			viewTarget.y = (int)(to.getY()+to.getHeight()-viewSize.getHeight()+
					SCROLL_MARGIN);
		}else{
			viewTarget.y = (int)(to.getY()-SCROLL_MARGIN);
		}		
		
		timer = new Timer(50,this);
		timer.start();
	}
	
	public void actionPerformed(ActionEvent e) {
		
		//If the previous execution of this method havent finished yet, just skip the frame
		if(!isBusy){
			isBusy = true;
			if(frame <= FRAME_COUNT){
				Point point = new Point();
				double factor = getFactor(frame);
				
				point.x = 
					(int)(viewTarget.getX() - (viewTarget.getX() - from.getX())*factor);
				point.y = 
					(int)(viewTarget.getY() - (viewTarget.getY() - from.getY())*factor);

				scrollable.setViewPosition(point);
			} else {
				//Make sure that we have reached the target
				scrollable.setViewPosition(viewTarget);
				timer.stop();
				latestAnimator = null;
			}
			isBusy = false;
		}
		frame++;
	}
	
	private void kill(){
		if(timer != null){
			timer.stop();
		}
		latestAnimator = null;
		logger.debug("ScrollAnimator killed");
	}
	
	/**
	 * This decides at what speed the scrolling is done. getFactor(1) should be about 1
	 * and getFactor(FRAME_COUNT) exactly 0. This one is slower at the end, which
	 * looks nice.
	 * 
	 * @param frame
	 * @return
	 */
	private double getFactor(double frame){
		return 1.0/frame - 1.0/FRAME_COUNT;
	}
}