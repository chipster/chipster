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
	
	private static final int COLOR_MARKER_WIDTH = 8;
	//To allow changing colors for the text
	private Color[] textColors;
	private String[] texts;
	private boolean emphasis = false;
    
    
    /**    
     * Text with multiple colors
     * 
     * @param x
     * @param y
     * @param z
     * @param texts Array of text which is shown concatenated
     * @param textColors Array of colors, indexes correspond to the indexes in array texts.
     */
    public Text(double x, double y, double z, String[] texts, Color[] textColors){
    	
        this.visualisationCoords = new double[1][3];
        this.projectedCoords = new double[1][2];
        this.deviceCoords = new int[1][2];
        this.dataCoords = new double[1][3];
        this.distanceFromCamera = new double[1];
    	
    	dataCoords[0][0] = visualisationCoords[0][0] = x;
        dataCoords[0][1] = visualisationCoords[0][1] = y;
        dataCoords[0][2] = visualisationCoords[0][2] = z;
    	
    	this.texts = texts;
    	this.textColors = textColors;
    }
    
    public Text(double x, double y, double z, String text, boolean emphasis) {
		this(x, y, z, new String[] {text}, new Color[] {Color.gray});
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
        
        if (this.texts != null) {
        	
        	//Calculate text dimensions
        	int textWidth = 0;
        	int textHeight = 10;

        	for(int i = 0; i < texts.length && i < textColors.length ; i++){
        		if (texts[i] != null) {
        			text = texts[i];        		
        			textWidth += SwingUtilities.computeStringWidth(g.getFontMetrics(), text);
        			if (!Color.gray.equals(textColors[i])) {
        				textWidth += COLOR_MARKER_WIDTH;
        			}
        		}
        	}
        	
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

        	int xShift = 0;        	        	
        	for(int i = 0; i < texts.length && i < textColors.length ; i++){
        		if (texts[i] != null) {
        			int textPartWidth = 0;
        			if (!Color.gray.equals(textColors[i])) {
        				g.setColor(textColors[i]);
        				g.fillRect(x + xShift, y - 6, 6, 6);
        				textPartWidth += COLOR_MARKER_WIDTH;
        			}
        			g.setColor(Color.gray);
        			Font origFont = g.getFont();
        			if (emphasis) {
        				g.setFont(g.getFont().deriveFont(Font.BOLD, 16));
        			}
        			g.drawString(texts[i], x + xShift + textPartWidth, y);
        			g.setFont(origFont);
        			
        			textPartWidth += SwingUtilities.computeStringWidth(g.getFontMetrics(), texts[i]);
        			
        			xShift += textPartWidth;
        		}
        	}
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
