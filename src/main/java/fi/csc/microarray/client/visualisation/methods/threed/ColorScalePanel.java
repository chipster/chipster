package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.JPanel;

/**
 * JPanel with black background and rainbow gradient. Values of the different colors
 * are taken from the dataModel and their colors from the color model.
 * 
 * @author Petri Klemelä
 */
public class ColorScalePanel extends JPanel{
	
	DataModel dataModel;
	float[] scaleValues;
	
	public ColorScalePanel(DataModel dataModel){
		this.dataModel = dataModel;
		this.scaleValues = dataModel.getColorScaleValues();
		this.setSize(50, 300);
		this.setPreferredSize(new Dimension(50, 300));
		this.setBackground(Color.BLACK);
		this.setOpaque(true);		
	}
	
	@Override
	protected void paintComponent(Graphics g){
		super.paintComponent(g);
		float upperLimit = (float)g.getClipBounds().getHeight() - 20;
		for(int y = 0 ; y < upperLimit ; y++){
			g.setColor(dataModel.getColorModel().getColorFor((y/upperLimit)));
			g.drawLine(10, 10+y, 20, 10+y);
		}
		
		for(float value : scaleValues){
			int y = (int)(((value - scaleValues[0]) / 
					(scaleValues[scaleValues.length-1] - scaleValues[0])) * upperLimit-1);
			
			//Little fix for rounding problem, which causes first line to be one pixel too high
			if(value == scaleValues[0]){
				y = 0;
			}
			
			g.setColor(Color.WHITE);
			
			int[] xPoints = new int[]{
					21, 26, 26
			};
			
			int[] yPoints = new int[]{
					10+y, 7+y, 13+y
			};
			
			g.fillPolygon(xPoints, yPoints, 3);
			//g.drawLine(21, 10+y, 25, 10+y);
			g.drawString(""+value, 30, y+15);
		}
	}
}
