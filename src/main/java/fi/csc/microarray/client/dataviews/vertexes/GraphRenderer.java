package fi.csc.microarray.client.dataviews.vertexes;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.Stroke;

import javax.swing.ImageIcon;

import org.jgraph.JGraph;
import org.jgraph.graph.CellView;
import org.jgraph.graph.VertexRenderer;

import fi.csc.microarray.client.VisualConstants;

/**
 * Renderer for a vertexes in MicroarrayGraph. For a collapsed group cell 
 * this renderer paints the small '+' and '-' handles in an upper left corner. 
 * It also takes care of rendering special borders if cell is selected.
 * 
 * @author mkoski
 */
public class GraphRenderer extends VertexRenderer {
	
	public static Rectangle phenodata = new Rectangle(0,0,0,0);
	
	/**
	 * Specifies whether the current view is a rich text value, and if the image
	 * should be stretched.
	 */
	protected boolean isGroup = false;
	
	protected boolean isPhenodata = false;

	protected boolean isSelected = false;


	protected Color graphForeground = Color.black;

	private JGraph graph;
	
	private AbstractGraphVertex vertex;

	/**
	 * Overrides the parent implementation to return the value component stored
	 * in the user object instead of this renderer if a value component exists.
	 * This applies some of the values installed to this renderer to the value
	 * component (border, opaque) if the latter is a JComponent.
	 * 
	 * @return Returns a configured renderer for the specified view.
	 */
	@Override
	public Component getRendererComponent(JGraph graph, CellView view,
			boolean sel, boolean focus, boolean preview) {
		isSelected = sel;
		graphForeground = graph.getForeground();
		isGroup = view.getCell() instanceof GroupVertex;
		isPhenodata = view.getCell() instanceof PhenodataVertex;
		if(view.getCell() instanceof AbstractGraphVertex){
			if(isGroup){
				vertex = (GroupVertex)view.getCell();
			} else if(isPhenodata){
				vertex = (PhenodataVertex)view.getCell();
			} else {
				vertex = (GraphVertex)view.getCell();
			}
		}
		this.graph = graph;
		
		return super.getRendererComponent(graph, view, sel, focus, preview);
	}

	/**
	 * Renderer paint method. Works like a default renderer but adds some 
	 * special effects it cell is selected. Paints also the groups 
	 * folding handle
	 */
	public void paint(Graphics g) {
		// Sets the font depending if the cell is selected
		if (isGroup) {
			if (hasFocus || selected) {
				 this.setFont(this.getFont().deriveFont(Font.BOLD));
			}
		} else {
			if (selected) {
				 this.setFont(this.getFont().deriveFont(Font.BOLD));
			}
		}

		if(isPhenodata){
			// Draw an oval
			
			int b = borderWidth;
			Graphics2D g2 = (Graphics2D) g;
			Dimension d = getSize();
			if (super.isOpaque()) {
				g.setColor(super.getBackground());
				if (gradientColor != null && !preview) {
					setOpaque(false);
					g2.setPaint(new GradientPaint(0, 0, getBackground(),
							getWidth(), getHeight(), gradientColor, true));
				}
				g.fillOval(b - 1, b - 1, d.width - b, d.height - b);
			}
			
			setBorder(null);
			super.paint(g);
			
			// Draw small oval borders
			if(!selected){
				g.setColor(bordercolor);
				g2.setStroke(new BasicStroke(b));
				g.drawOval(0, 0, d.width - 1, d.height - 1);
			}
		} else {
			super.paint(g);
		}		
		
		if(isPhenodata && !((PhenodataVertex)vertex).isPhenodataSet()){
			// Draws warning icon for phenodata vertex if phenodata is not set
			ImageIcon icon = VisualConstants.PHENODATA_ICON;
			
			double iconHeight = icon.getIconHeight();
			double iconWidth = icon.getIconWidth();
			Dimension d = getSize();
			Image img = icon.getImage();
		
			// Icons position is calculated using these factors
			double xPos = 2.4f;
			double yPos = 1.3f;
			
			// Increase the clip area. This must be done, because the icon is 
			// drawn outside the original clip area.
			g.setClip(
					(int)(g.getClipBounds().getX() - iconWidth / xPos), 				// x
					(int)(g.getClipBounds().getY()), 									// y
					(int)(g.getClipBounds().getWidth() + iconWidth / (1/xPos)), 		// width
					(int)(g.getClipBounds().getHeight() + iconHeight / (1/yPos))		// height
					);
			g.drawImage(
					img, 
					(int)(-iconWidth / xPos), 
					(int)(d.getHeight() - iconHeight / yPos), 
					(int)iconWidth, 
					(int)iconHeight, 
					graph);
			
			phenodata.setBounds(new Rectangle(					
					(int)(-iconWidth / xPos), 
					(int)(d.getHeight() - iconHeight / yPos), 
					(int)iconWidth, 
					(int)iconHeight));
		}
	}

	/**
	 * Paint special borders for selected cells. 
	 * This method is called by the paint method.
	 */
	@Override
	public void paintSelectionBorder(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		Stroke previousStroke = g2.getStroke();

		boolean paintBorders = false;
			// Border for a normal vertex which is selected
			if (selected) {
				g2.setStroke(GraphVertex.SELECTION_STROKE);
				g.setColor(GraphVertex.SELECTED_BORDER_COLOR);
				paintBorders = true;
			}		

		// Paints the border if needed
		if (paintBorders) {
			Dimension d = getSize();
			if(isPhenodata){
				// Draw oval borders
				g.drawOval(0, 0, d.width - 1, d.height - 1);
			} else {
				g.drawRect(0, 0, d.width - 1, d.height - 1);
			}
		}

		// Sets the stroke as if was before this method
		g2.setStroke(previousStroke);
	}
}
