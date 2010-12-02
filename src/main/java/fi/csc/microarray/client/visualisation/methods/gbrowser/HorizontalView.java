package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.geom.Point2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegionDouble;

/**
 * The basic view of genome browser. 
 *
 */
public class HorizontalView extends View implements KeyListener {

	public HorizontalView(GenomePlot parent, boolean movable, boolean zoomable, boolean selectable) {
				
		super(parent, movable, zoomable, selectable);
		
		parent.chartPanel.addKeyListener(this);
	}

	@Override
	protected void drawView(Graphics2D g, boolean isAnimation) {

		super.drawView(g, isAnimation);

		// show current position on top of chromosome
		if (highlight != null) {
			g.setPaint(new Color(0, 0, 0, 64));
			Rectangle rect = g.getClip().getBounds();

			rect.x = bpToTrack(highlight.start);
			rect.width = Math.max(1, bpToTrack(highlight.end) - rect.x);
			rect.height = 24;
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

		g.setFont(g.getFont().deriveFont(10f));
		TextDrawable text = (TextDrawable) drawable;

		text.text = text.text.replaceAll("\"", "");

		g.drawString(text.text, text.x + x, text.y + y);
	}

	protected void drawRectDrawable(Graphics2D g, int x, int y, Drawable drawable) {

		RectDrawable rect = (RectDrawable) drawable;

		// Draw fill
		if (rect.color != null) {
			g.setPaint(rect.color);
			g.fillRect(rect.x + x, rect.y + y, rect.width, rect.height);
		}

		// Draw outline after fill to make sure that it stays continuous
		if (rect.lineColor != null) {
			g.setPaint(rect.lineColor);
			g.drawRect(rect.x + x, rect.y + y, rect.width-1, rect.height-1);
		}
	}

	protected void drawLineDrawable(Graphics2D g, int x, int y, Drawable drawable) {
		LineDrawable line = (LineDrawable) drawable;
		g.drawLine(line.x + x, line.y + y, line.x2 + x, line.y2 + y);
	}

	@Override
	protected void handleDrag(Point2D start, Point2D end, boolean disableDrawing) {
		double bpMove = trackToBp((double) start.getX()).minus(trackToBp((double) end.getX()));

		bpRegion.move(bpMove);
		setBpRegion(bpRegion, disableDrawing);

		parentPlot.redraw();
	}
	
	@Override
	public void setBpRegion(BpCoordRegionDouble region, boolean disableDrawing) {
		
		BpCoord maxBp = getMaxBp();
		
		if (maxBp != null && region.getLength() > maxBp.bp) {
			region.start.bp = 0.0;
			region.end.bp = (double)maxBp.bp;
		} else {
			if (region.start.bp < 0 ) {
				region.move(-region.start.bp);
			}

			if (maxBp != null && region.end.bp > maxBp.bp) {				
				double delta = region.end.bp - (double)maxBp.bp;
				region.move(-delta);
			}
		}
		
		super.setBpRegion(region, disableDrawing);
	}

	@Override
	public void keyPressed(KeyEvent e) {
		
		if ( this.zoomable ) {
			if ( e.getKeyCode() == KeyEvent.VK_UP ) {
				
				this.zoomAnimation(this.getWidth() / 2, -1);
				
			} else if (e.getKeyCode() == KeyEvent.VK_DOWN )  {
				
				this.zoomAnimation(this.getWidth() / 2, 1);
			}
		}
		
		if ( this.movable ) {
			if ( e.getKeyCode() == KeyEvent.VK_LEFT ) {

				bpRegion.move(-this.getBpRegion().getLength() / 10.0);
				setBpRegion(bpRegion, false);
				parentPlot.redraw();
				
			} else if (e.getKeyCode() == KeyEvent.VK_RIGHT )  {

				bpRegion.move(this.getBpRegion().getLength() / 10.0);
				setBpRegion(bpRegion, false);
				parentPlot.redraw();
			}
		}
	}

	@Override
	public void keyReleased(KeyEvent e) {

	}

	@Override
	public void keyTyped(KeyEvent e) {
					
	}
}
