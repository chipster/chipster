package fi.csc.microarray.client.dataview;

import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Timer;

import fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex;

public class PositionChangeAnimator implements ActionListener{
	private Point from;
	private Point to;
	private int frame = 1;
	private AbstractGraphVertex cell;
	private Timer timer;
	private static boolean isBusy = false;

	public static final int FRAME_COUNT = 5;

	public PositionChangeAnimator(AbstractGraphVertex cell, Point from, Point to){
		this.cell = cell;
		this.from = from;
		this.to = to;
		
		timer = new Timer(50,this);
		timer.start();
	}
	public void actionPerformed(ActionEvent e) {
//		If the previous execution of this method havent finished yet, just skip the frame
		if(!isBusy){
			isBusy = true;
			if(frame <= FRAME_COUNT){
				Point point = new Point();
				double factor = getFactor(frame);
				point.x = 
					(int)(to.getX() - (to.getX() - from.getX())*factor);
				point.y = 
					(int)(to.getY() - (to.getY() - from.getY())*factor);

				cell.setPosition(point);
			} else {
				cell.setPosition(to);
				timer.stop();
			}
			isBusy = false;
		}
		frame++;
	}
	
	private double getFactor(double frame){
		return 1.0/frame - 1.0/FRAME_COUNT;
	}
}