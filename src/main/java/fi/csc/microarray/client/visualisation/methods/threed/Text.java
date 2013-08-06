package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;

/**
 * @author klemela
 */
public class Text extends Drawable {
	
	//To allow changing colors for the text
	private Color textColor;
	private String text;
	private boolean emphasis = false;
        
    /**    
     * Text with color
     * 
     * @param x
     * @param y
     * @param z
     * @param texts Array of text which is shown concatenated
     * @param textColors Array of colors, indexes correspond to the indexes in array texts.
     */
    public Text(float x, float y, float z, String text, Color color){
    	
        this.visualisationCoords = new double[1][3];
        this.projectedCoords = new double[1][2];
        this.deviceCoords = new int[1][2];
        this.dataCoords = new double[1][3];
        this.distanceFromCamera = new double[1];
    	
    	dataCoords[0][0] = visualisationCoords[0][0] = x;
        dataCoords[0][1] = visualisationCoords[0][1] = y;
        dataCoords[0][2] = visualisationCoords[0][2] = z;
    	
    	this.text = text;
    	this.textColor = color;
    }
    
    public Text(float x, float y, float z, String text, Color color, boolean emphasis) {
		this(x, y, z, text, color);
		this.emphasis = emphasis;
	}

	/**
     * 
     * @param g 
     * @param width 
     * @param height 
     */
    public void draw(Graphics g, int width, int height, PaintMode notUsed) {    
        
        g.setColor(color);
        deviceCoords[0][0] = (int)((projectedCoords[0][0] + 0.5) * width);
        deviceCoords[0][1] = (int)((projectedCoords[0][1] + 0.5) * height);
        
        if (this.text != null) {
        	
        	//Calculate text dimensions
        	int textWidth = 0;
        	int textHeight = 10;
   		
   			textWidth += SwingUtilities.computeStringWidth(g.getFontMetrics(), text);
        	
        	int textCenterX = deviceCoords[0][0];
        	int textCenterY = deviceCoords[0][1];
        	
        	int x;
        	int y;
        	
        	if (textCenterX < width / 2) {
        		x = textCenterX - textWidth;        		
        	} else {
        		x = textCenterX;
        	}
        	        	
        	if (textCenterY < height / 2) {        	
        		y = textCenterY;        
        	} else {
        		y = textCenterY + textHeight;
        	}

			g.setColor(textColor);
			
			Font origFont = g.getFont();
			if (emphasis) {
				g.setFont(g.getFont().deriveFont(Font.BOLD, 16));
			}
			g.drawString(text, x, y);
			g.setFont(origFont);        			        	
        }
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
    }
}
