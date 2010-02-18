package fi.csc.microarray.client.visualisation.methods.genomeBrowser;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Point2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;

public class VerticalView  extends View{

	public VerticalView(GenomeBrowser parent, boolean movable, boolean zoomable,
			boolean selectable) {
		super(parent, movable, zoomable, selectable);		
	}

	@Override
	protected void drawView(Graphics2D g, boolean isAnimation){

		super.drawView(g, isAnimation);
		
		if(highlight != null){

			g.setPaint(new Color(0,0,0,64));
			Rectangle rect = g.getClip().getBounds();
			rect.y = bpToTrack(highlight.start);
			rect.height = bpToTrack(highlight.end) - rect.y;
			g.fill(rect);
		}
	}

	protected void drawDrawable(Graphics2D g, int x, int y, Drawable drawable) {
		
		int tmp = x;
		y = x;
		x = tmp;

		g.setPaint(drawable.color);			
		// TODO Auto-generated method stub

		if (drawable instanceof TextDrawable) {

			TextDrawable text = (TextDrawable) drawable;
			g.drawString(text.text, text.y + y + 10, text.x + x);		

		} else if (drawable instanceof RectDrawable) {

			RectDrawable rect = (RectDrawable) drawable;

			if(rect.lineColor != null){
				g.setPaint(rect.lineColor);			
				g.fillRect(rect.y + y, rect.x + x, rect.height, rect.width);
			}
			//Draw outline after fill to make sure that it stays continuous
			g.setPaint(drawable.color);
			g.drawRect(rect.y + y, rect.x + x, rect.height, rect.width);

		} else if (drawable instanceof LineDrawable) {

			LineDrawable line = (LineDrawable) drawable;
			g.drawLine(line.y + y, line.x + x, line.y2 + y, line.x2 + x);
		}
	}

	protected void handleDrag(Point2D start, Point2D end, boolean disableDrawing) {
		double bpMove = trackToBp((double)start.getY()).minus(trackToBp((double)end.getY()));

		if(bpMove < 0 && bpRegion.start.bp < Math.abs(bpMove)){
			bpMove = -bpRegion.start.bp;
		}

		bpRegion.move(bpMove);
		setBpRegion(bpRegion, disableDrawing);		

		parentPlot.redraw();
	}
	
	@Override
	public float getTrackHeight() {
		 return this.getHeight() / (float)tracks.size();
	}

	@Override
	public int getWidth() {
		return super.getHeight();
	}

	@Override
	public int getHeight() {
		return super.getWidth();
	}
	
	public int getX() {
		return super.getY();
	}
}
