package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Point2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;

/**
 * The basic view of genome browser. 
 *
 */
public class HorizontalView extends View {

	public HorizontalView(GenomePlot parent, boolean movable, boolean zoomable, boolean selectable) {
		super(parent, movable, zoomable, selectable);
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

		if (bpMove < 0 && bpRegion.start.bp < Math.abs(bpMove)) {
			bpMove = -bpRegion.start.bp;
		}
		
		BpCoord maxBp = getMaxBp();
		
		if (maxBp != null && bpRegion.end.bp + bpMove > maxBp.bp) {
			bpMove = maxBp.bp - bpRegion.end.bp;
		}

		bpRegion.move(bpMove);
		setBpRegion(bpRegion, disableDrawing);

		parentPlot.redraw();
	}
}
