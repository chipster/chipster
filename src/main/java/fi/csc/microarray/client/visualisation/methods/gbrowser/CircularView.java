package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.MouseWheelEvent;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

public class CircularView extends View {

	private Point offset;
	private float angleOffset = (float) Math.PI;

	public CircularView(GenomePlot parent, boolean movable, boolean zoomable, boolean selectable) {
		super(parent, movable, zoomable, selectable);
	}

	@Override
	protected void drawView(Graphics2D g, boolean isAnimation) {

		super.drawView(g, isAnimation);

		if (highlight != null) {

			g.setPaint(new Color(0, 0, 0, 64));
			Rectangle rect = g.getClip().getBounds();
			rect.x = bpToTrack(highlight.start);
			rect.width = bpToTrack(highlight.end) - rect.x;
			g.fill(rect);
		}
	}

	@Override
	protected void drawDrawable(Graphics2D g, int x, int y, Drawable drawable) {

		this.offset = new Point(x, y);

		g.setPaint(drawable.color);

		if (drawable instanceof TextDrawable) {

			TextDrawable text = (TextDrawable) drawable;

			int pixelLength = 0;

			for (int i = 0; i < text.text.length(); i++) {

				char chr = text.text.charAt(i);
				Point.Float p = toCircle(new Point(text.x + pixelLength, text.y - 4));

				if (p != null) {

					TextLayout tl = new TextLayout("" + text.text.charAt(i), g.getFont(), g.getFontRenderContext());

					Shape rotatedChar = tl.getOutline(AffineTransform.getRotateInstance(getAngle(p), tl.getBounds().getCenterX(), tl.getBounds().getCenterY()));

					g.translate(p.x - tl.getBounds().getCenterX(), p.y - tl.getBounds().getCenterY());
					g.fill(rotatedChar);
					g.translate(-(p.x - tl.getBounds().getCenterX()), -(p.y - tl.getBounds().getCenterY()));
				}

				pixelLength += g.getFontMetrics().charWidth(chr);
			}

		} else if (drawable instanceof RectDrawable) {

			RectDrawable rect = (RectDrawable) drawable;

			// Fill
			drawPerimeterLine(g, rect);

			// Radial outlines
			drawRectRadialOutline(g, rect);

			// Outlines along perimeter
			int height = rect.height;
			rect.height = 1;
			rect.color = rect.lineColor;
			drawPerimeterLine(g, rect);
			rect.y += height;
			drawPerimeterLine(g, rect);

		} else if (drawable instanceof LineDrawable) {

			LineDrawable line = (LineDrawable) drawable;
			Point.Float p1 = new Point.Float(line.x, line.y);
			Point.Float p2 = new Point.Float(line.x2, line.y2);
			drawAnyLine(g, p1.x, p1.y, p2.x, p2.y);
		}
	}

	private float asDeg(float radian) {
		return (float) (radian * 360 / (Math.PI * 2));
	}

	private void drawRectRadialOutline(Graphics2D g, RectDrawable rect) {

		g.setPaint(rect.lineColor);

		float x1 = rect.x;
		float x2 = rect.x + rect.width;
		float y1 = rect.y;
		float y2 = rect.y + rect.height;

		drawAnyLine(g, x1, y1, x1, y2);
		drawAnyLine(g, x2, y1, x2, y2);
	}

	private void drawPerimeterLine(Graphics2D g, RectDrawable rect) {

		rect.height = Math.max(rect.height, 0);

		g.setPaint(rect.color);
		g.setStroke(new BasicStroke(rect.height, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));

		Point.Float startPoint = toCircle(new Point.Float(rect.x, rect.y));
		Point.Float endPoint = toCircle(new Point.Float(rect.x + rect.width, rect.y));

		// Null means that angle is over Math.PI * 2, this drawable should be clipped
		if (startPoint != null && endPoint != null) {

			float r = rect.y + offset.y + rect.height / 2;

			float startAngle = getAngle(startPoint);

			float start = asDeg(getAngle(startPoint));
			float extent = asDeg(getAngle(endPoint) - startAngle);

			// Angles of Arc2D start from three o'clock and grow counter-clockwise
			start = -(start - 90);
			extent = -extent;

			g.draw(new Arc2D.Float(getCenterX() - r, getCenterY() - r, r * 2, r * 2, start, extent, Arc2D.OPEN));
		}
		g.setStroke(new BasicStroke(1));

	}

	private void drawAnyLine(Graphics2D g, float x1, float y1, float x2, float y2) {

		// Long lines are drawn in parts to make them curve, +1 to make sure there is at least one
		// part
		int parts = (int) (Math.abs(x2 - x1) * (Math.abs(getAngle(x2) - getAngle(x1))) + 1);

		float dx = (x2 - x1) / parts;
		float dy = (y2 - y1) / parts;

		Point.Float p2 = null;

		for (int i = 0; i <= parts; i++) {
			Point.Float p1 = toCircle(new Point.Float(x1, y1 + 0.5f));

			if (p2 != null && p1 != null) { // p2 null on first round, p1 null if goes over clip limits
				g.drawLine((int) p1.x, (int) p1.y, (int) p2.x, (int) p2.y);
			}
			p2 = p1;
			x1 += dx;
			y1 += dy;
		}
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		if (zoomable) {

			double pointerBp = pointToBp(new Point.Float(e.getPoint().x, e.getPoint().y)).bp;
			double pointerRelative = (pointerBp - bpRegion.start.bp) / bpRegion.getLength();

			if (e.getWheelRotation() > 0) {
				pointerRelative = 1 - pointerRelative;
			}

			pointerBp = bpRegion.start.bp + pointerRelative * bpRegion.getLength();

			double startBp = getBpRegionDouble().start.bp;
			double endBp = getBpRegionDouble().end.bp;

			double width = endBp - startBp;
			width *= Math.pow(ZOOM_FACTOR, e.getWheelRotation());

			startBp = (double) (pointerBp - width * pointerRelative);
			endBp = (double) (pointerBp + width * (1 - pointerRelative));

			if (startBp < 0) {
				endBp += -startBp;
				startBp = 0;
			}

			setBpRegion(new BpCoordRegionDouble(startBp, endBp, new Chromosome("1")), false);
		}
	}

	protected void handleDrag(Point2D start, Point2D end, boolean disableDrawing) {

		Point.Float startPoint = new Point.Float((float) start.getX(), (float) start.getY());
		Point.Float endPoint = new Point.Float((float) end.getX(), (float) end.getY());

		float a1 = getAngle(startPoint);
		float a2 = getAngle(endPoint);

		double bpMove = pointToBp(startPoint).minus(pointToBp(endPoint));

		if (a1 > Math.PI * (3 / 2) && a2 < Math.PI / 2) {
			bpMove -= bpRegion.getLength();
		}

		if (a2 > Math.PI * (3 / 2) && a1 < Math.PI / 2) {

			bpMove += bpRegion.getLength();
		}

		if (bpMove < 0 && bpRegion.start.bp < Math.abs(bpMove)) {
			bpMove = -bpRegion.start.bp;
		}

		bpRegion.move(bpMove);
		setBpRegion(bpRegion, disableDrawing);

		parentPlot.redraw();
	}

	private Point.Float toCircle(Point xy) {
		return toCircle(new Point.Float(xy.x, xy.y));
	}

	private Point.Float toCircle(Point.Float xy) {

		float r = xy.y + offset.y;

		if (xy.x < offset.x || xy.x > this.getWidth() + offset.x - margin) {
			return null;
		}

		float angle = -getAngle(xy.x - offset.x);
		float dx = (int) (Math.sin(angle) * r);
		float dy = (int) (Math.cos(angle) * r);

		return new Point.Float(dx + getCenterX(), dy + getCenterY());
	}

	private int getCenterX() {
		return getWidth() / 2 + offset.x;
	}

	private int getCenterY() {
		return getHeight();
	}

	@Override
	public int getHeight() {
		return super.getHeight() / 2;
	}

	private float getAngle(Point.Float p) {

		float a = (float) (Math.atan((p.y - getCenterY()) / (p.x - getCenterX())) + Math.PI / 2);

		if (p.x < getCenterX()) {
			a += Math.PI;
		}
		return a;
	}

	private BpCoordDouble pointToBp(Point.Float p) {
		return trackToBp((double) (getAngle(p) * getWidth() / 2 / Math.PI + offset.x));
	}

	private float getAngle(float x) {
		return x / (float) getWidth() * (float) Math.PI * 2 + angleOffset;
	}

	// @Override
	// public float getTrackHeight() {
	// if(trackHeight == null){
	// trackHeight = (getHeight() - getMaxTrackHeightTotal() * 2) / ((float)getStretchableTrackCount() * 2);
	// }
	// return trackHeight;
	// }
}
