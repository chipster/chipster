package fi.csc.microarray.client.visualisation.methods.genomeBrowser;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Point2D;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.TextDrawable;

public class HorizontalView  extends View{

	public HorizontalView(GenomeBrowser parent, boolean movable,
			boolean zoomable, boolean selectable) {
		super(parent, movable, zoomable, selectable);
	}

	@Override
	protected void drawView(Graphics2D g, boolean isAnimation){
		
		super.drawView(g, isAnimation);

		if(highlight != null){

			g.setPaint(new Color(0,0,0,64));
			Rectangle rect = g.getClip().getBounds();
			rect.x = bpToTrack(highlight.start);
			rect.width = Math.max(1, bpToTrack(highlight.end) - rect.x);
			g.fill(rect);
		}
	}

	protected void drawDrawable(Graphics2D g, int x, int y, Drawable drawable) {

		g.setPaint(drawable.color);			

		if (drawable instanceof TextDrawable) {
			drawTextDrawable(g, x, y, drawable);

		} else if (drawable instanceof RectDrawable) {
			drawRectDrawable(g, x, y, drawable);

		} else if (drawable instanceof LineDrawable) {
			drawLineDrawable(g, x, y, drawable);
		}
	}
	
	protected void drawTextDrawable(Graphics2D g, int x, int y, Drawable drawable) {
		
		TextDrawable text = (TextDrawable) drawable;
		g.drawString(text.text, text.x + x, text.y + y);		
	}
	
	protected void drawRectDrawable(Graphics2D g, int x, int y, Drawable drawable) {
		
		RectDrawable rect = (RectDrawable) drawable;
		
		if(rect.color != null){
			
			if(rect.lineColor == null){
				rect.x -= 1;
				//rect.y -= 1;
				rect.width += 2;
				//rect.height += 2;
			}
			
			g.setPaint(drawable.color);			
			g.fillRect(rect.x + x, rect.y + y, rect.width, rect.height);
			
		}
		
		//Draw outline after fill to make sure that it stays continuous			
		if(drawable.color != null){
			g.setPaint(rect.lineColor);
			g.drawRect(rect.x + x, rect.y + y, rect.width, rect.height);
		}
	}
	
	protected void drawLineDrawable(Graphics2D g, int x, int y, Drawable drawable) {
		LineDrawable line = (LineDrawable) drawable;
		g.drawLine(line.x + x, line.y + y, line.x2 + x, line.y2 + y);		
	}
	
	

	@Override
	protected void handleDrag(Point2D start, Point2D end, boolean disableDrawing) {
		long bpMove = trackToBp((long)start.getX()) - trackToBp((long)end.getX());

		if(bpMove < 0 && bpRegion.start < Math.abs(bpMove)){
			bpMove = -bpRegion.start;
		}

		bpRegion.move(bpMove);
		setBpRegion(bpRegion, disableDrawing);		

		parentPlot.redraw();
	}
}
