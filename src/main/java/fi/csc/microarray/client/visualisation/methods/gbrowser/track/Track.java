package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Point;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;

public abstract class Track implements AreaResultListener {

	protected View view;
	protected DataSource file;
	protected Strand strand = Strand.FORWARD;

	public Track(View view, DataSource file) {
		this.view = view;
		this.file = file;
	}

	public Track(View view, DataSource file, Class<? extends AreaRequestHandler> handler, FileParser inputParser) {
		this(view, file);
		view.getQueueManager().createQueue(file, handler, inputParser);
	}

	/**
	 * Should be called after Track object is created, but can't be merged to constructor, because the coming areaResult event could cause
	 * call to track object before it's constructed.
	 */
	public void initializeListener() {
		if (file != null) {
			view.getQueueManager().addResultListener(file, this);
		} else {
			System.out.println("Track has no file: " + this);
		}
	}

	public abstract Collection<Drawable> getDrawables();

	protected View getView() {
		return view;
	}

	public void updateData() {
		if (file != null && this.getMaxHeight() > 0) {
			FsfStatus status = new FsfStatus();
			status.clearQueues = true;
			status.concise = isConcised();

			Collection<ColumnType> defCont = getDefaultContents();
			
			BpCoordRegion requestRegion = new BpCoordRegion(view.getBpRegion());
			requestRegion.start.bp -= 1000;
			requestRegion.end.bp += 1000;
			
			if (requestRegion.start.bp < 0) {
				requestRegion.start.bp = 0l;
			}

			view.getQueueManager().addAreaRequest(file, new AreaRequest(requestRegion, defCont, status), true);
		}
	}

	public abstract Collection<ColumnType> getDefaultContents();

	public abstract boolean isConcised();

	public Collection<Drawable> getEmptyDrawCollection() {
		return new LinkedList<Drawable>();
	}

	public int getMaxHeight() {
		return Integer.MAX_VALUE;
	}

	public void setStrand(Strand s) {
		this.strand = s;
	}

	public Strand getStrand() {
		return strand;
	}

	private Point2D[] arrowPoints = new Point2D[] { new Point.Double(0, 0.25), new Point.Double(0.5, 0.25), new Point.Double(0.5, 0), new Point.Double(1, 0.5), new Point.Double(0.5, 1), new Point.Double(0.5, 0.75), new Point.Double(0, 0.75), new Point.Double(0, 0.25) };

	protected Collection<? extends Drawable> getArrowDrawables(int x, int y, int width, int height) {

		Collection<Drawable> parts = getEmptyDrawCollection();

		for (int i = 1; i < arrowPoints.length; i++) {
			Point2D p1 = arrowPoints[i - 1];
			Point2D p2 = arrowPoints[i];

			Point2D p1Scaled = new Point.Double(x + p1.getX() * width, y + p1.getY() * height);
			Point2D p2Scaled = new Point.Double(x + p2.getX() * width, y + p2.getY() * height);

			parts.add(new LineDrawable((int) p1Scaled.getX(), (int) p1Scaled.getY(), (int) p2Scaled.getX(), (int) p2Scaled.getY(), Color.black));
		}

		return parts;
	}

	public BpCoord getMaxBp(Chromosome chr) {
		return null;
	}
}
