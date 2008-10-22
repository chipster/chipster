package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.Icon;

/**
 * 
 * Simple class implementing Icon interface, which dwaws colored circles. This
 * is used in the OperationCategory list to show different color for every
 * category.
 * 
 * @author klemela
 *
 */
public class ColoredCircleIcon implements Icon {
	private int height, width;
	private Color color;
	
	public ColoredCircleIcon(Color c){
		color = c;
		this.height = 8;
		this.width = 8;
	}
	
	public ColoredCircleIcon(Color c, int diameter){
		color = c;
		this.height = diameter;
		this.width = diameter;
	}

	public int getIconHeight() {
		return height;
	}

	public int getIconWidth() {
		return width;
	}

	public void paintIcon(Component comp, Graphics g, int x, int y) {
		Graphics2D g2d = (Graphics2D)g;
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setColor(color);
	    g2d.fillOval(x,y, height, width);
	    //This draws one pixel too pig circle which causes strange effects on different environments (windows)
	    //g2d.drawOval(x,y, height, width);
	}

}
