	package fi.csc.microarray.client.visualisation.methods.gbrowser;

	import java.awt.Color;
	import java.awt.Font;
	import java.awt.Rectangle;
	import java.awt.Shape;
	import java.io.FileNotFoundException;
	import java.io.IOException;
	import java.io.ObjectInputStream;
	import java.io.ObjectOutputStream;
	import java.io.Serializable;
	import java.net.MalformedURLException;
	import java.util.Collection;
	import java.util.LinkedList;
	import java.util.List;

	import org.jfree.chart.ChartMouseEvent;
	import org.jfree.chart.ChartMouseListener;
	import org.jfree.chart.ChartPanel;
	import org.jfree.chart.plot.Plot;
	import org.jfree.chart.plot.PlotRenderingInfo;
	import org.jfree.chart.plot.PlotState;
	import org.jfree.data.general.DatasetChangeEvent;
	import org.jfree.util.ObjectUtilities;

	import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
	import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegionDouble;
	import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
	import fi.csc.microarray.client.visualisation.methods.gbrowser.track.EmptyTrack;

	/**
	 * @author Petri Klemelä, Aleksi Kallio
	 */
	public class GenomePlot extends Plot implements ChartMouseListener, Cloneable, Serializable { // , MouseWheelListener {

		//private static final File FILE_ROOT = new File("/home/akallio/chipster-share/genome_browser");
		private static final File FILE_ROOT = new File("/home/klemela/chipster-share/genome_browser");

//		private static final URL URL_ROOT;
	//
//		static {
//			try {
//				URL_ROOT = new URL("http://chipster-devel.csc.fi:8050/public/space_separated_annotations");
//			} catch (MalformedURLException e) {
//				throw new RuntimeException(e);
//			}
//		}

		/** The cell information text is drawn with this font. */
		private Font descriptionFont;

		private List<View> views = new LinkedList<View>();
		private View dataView = null;

		public View getDataView() {
			return dataView;
		}

		public View getOverviewView() {
			return overviewView;
		}

		private View overviewView = null;

		public ChartPanel chartPanel;

		public GenomePlot(boolean horizontal) throws FileNotFoundException, MalformedURLException {

			// add overview view
			this.overviewView = new HorizontalView(this, false, false, true);
			this.overviewView.margin = 0;
			this.views.add(overviewView);

			// add horizontal or circular data view
			if (horizontal) {
				this.dataView = new HorizontalView(this, true, true, false);

			} else {
				this.dataView = new CircularView(this, true, true, false);
				this.dataView.margin = 20;
				this.dataView.addTrack(new EmptyTrack(dataView, 30));
			}

			this.views.add(dataView);

			dataView.addRegionListener(new RegionListener() {
				public void RegionChanged(BpCoordRegion bpRegion) {
					overviewView.highlight = bpRegion;
				}
			});
		}

//		public void jee(boolean horizontal, boolean transcripts, boolean genes, boolean mirna, boolean sequence) throws FileNotFoundException, MalformedURLException {
	//
//			// cytobands
//			TrackFactory.addCytobandTracks(overviewView, dataView, new DataSource(URL_ROOT, "cytoband_hg17_sorted.fsf"));
	//
//			// Reference genes
//			if (genes) {
//				TrackFactory.addGeneTracks(this, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_genes.fsf"));
//			}
	//
//			// miRNA genes
//			if (mirna) {
//				TrackFactory.addMirnaTracks(dataView, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
//			}
	//
//			// Reference transcripts
//			if (transcripts) {
//				TrackFactory.addTranscriptTracks(dataView, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_transcripts.fsf"));
//			}
	//
//			// Peaks
//			TrackFactory.addPeakTracks(dataView, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
	//
//			// Wiggle
//			TrackFactory.addWigTrack(dataView, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
	//
//			// Eland export
//			TrackFactory.addReadTracks(dataView, new DataSource(FILE_ROOT, "eland_result.fsf"), new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_seq.fsf"));
	//
//			TrackFactory.addRulerTrack(dataView);
	//
//		}

		public void start(String chromosome) {
			dataView.setBpRegion(new BpCoordRegionDouble(0d, 1024 * 1024 * 250d, new Chromosome(chromosome)), false);
			overviewView.setBpRegion(new BpCoordRegionDouble(0d, 1024 * 1024 * 250d, new Chromosome(chromosome)), false);
		}
		public String getPlotType() {
			return "GeneBrowser";
		}

		public Font getDescriptionFont() {
			return this.descriptionFont;
		}

		/**
		 * Draws the plot on a Java2D graphics device (such as the screen or a printer).
		 * 
		 * @param g2
		 *            the graphics device.
		 * @param area
		 *            the area within which the plot should be drawn.
		 * @param anchor
		 *            the anchor point (<code>null</code> permitted).
		 * @param parentState
		 *            the state from the parent plot, if there is one (<code>null</code> permitted.)
		 * @param info
		 *            collects info about the drawing (<code>null</code> permitted).
		 * @throws NullPointerException
		 *             if g2 or area is null.
		 */
		public void draw(java.awt.Graphics2D g2, java.awt.geom.Rectangle2D area, java.awt.geom.Point2D anchor, PlotState parentState, PlotRenderingInfo info) {

			// adjust for insets...
			// this.getInsets().trim(area);

			if (info != null) {
				info.setPlotArea(area);
				info.setDataArea(area);
			}

			// this.setBackgroundPaint(Color.black);

			drawBackground(g2, area);
			drawOutline(g2, area);

			Shape savedClip = g2.getClip();
			g2.clip(area);
			/*
			 * Composite originalComposite = g2.getComposite(); g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,
			 * getForegroundAlpha()));
			 */

			Rectangle viewArea = (Rectangle) area.getBounds().clone();

			// Horizontal or vertical split
			if (true) {

				float[] viewHeights = new float[] { 0.1f, 0.9f };

				for (int i = 0; i < views.size(); i++) {
					if (i > 0) {
						viewArea.y += viewArea.height;
					}
					viewArea.height = (int) (area.getBounds().getHeight() * viewHeights[i]);
					g2.setClip(viewArea);
					views.get(i).drawView(g2, false);
				}

			} else {
				float[] viewWidths = new float[] { 0.05f, 0.95f };
				Rectangle lastArea = null;

				for (int i = 0; i < views.size(); i++) {
					if (lastArea != null) {
						viewArea.x = lastArea.x + lastArea.width;
						viewArea.y = lastArea.y;
						viewArea.height = lastArea.height;
					}
					g2.setColor(Color.black);
					viewArea.width = (int) (area.getBounds().getWidth() * viewWidths[i]);
					lastArea = (Rectangle) (viewArea.clone());

					View view = views.get(i);

					if (view instanceof VerticalView) {
						viewArea.grow(0, -view.margin);
					} else if (view instanceof HorizontalView) {
						viewArea.grow(-view.margin, 0);
					}

					g2.setClip(savedClip);
					g2.drawLine(viewArea.x - 1, 0, viewArea.x - 1, viewArea.height);

					g2.setClip(viewArea);
					view.drawView(g2, false);
				}
			}
			g2.setClip(savedClip);
			// g2.setComposite(originalComposite);

			drawOutline(g2, area);
		}

		/**
		 * Implements the ChartMouseListener interface. This method does nothing.
		 * 
		 * @param event
		 *            the mouse event.
		 */
		public void chartMouseMoved(ChartMouseEvent event) {
		}

		public void chartMouseClicked(ChartMouseEvent e) {
		}

		/**
		 * Tests this plot for equality with an arbitrary object. Note that the plot's dataset is NOT included in the test for equality.
		 * 
		 * @param obj
		 *            the object to test against (<code>null</code> permitted).
		 * 
		 * @return <code>true</code> or <code>false</code>.
		 */
		public boolean equals(Object obj) {
			if (obj == this) {
				return true;
			}
			if (!(obj instanceof GenomePlot)) {
				return false;
			}
			GenomePlot that = (GenomePlot) obj;

			// if (!ObjectUtilities.equal(this.urlGenerator, that.urlGenerator)) {
			// return false;
			// }

			if (!ObjectUtilities.equal(this.descriptionFont, that.descriptionFont)) {
				return false;
			}

			// can't find any difference...
			return true;
		}

		/**
		 * Provides serialization support.
		 * 
		 * @param stream
		 *            the output stream.
		 * 
		 * @throws IOException
		 *             if there is an I/O error.
		 * @throws NullPointerException
		 *             if stream is null.
		 */
		private void writeObject(ObjectOutputStream stream) throws IOException {
			stream.defaultWriteObject();
		}

		/**
		 * Provides serialization support.
		 * 
		 * @param stream
		 *            the input stream.
		 * 
		 * @throws IOException
		 *             if there is an I/O error.
		 * @throws ClassNotFoundException
		 *             if there is a classpath problem.
		 * @throws NullPointerException
		 *             if stream is null.
		 */
		private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
			stream.defaultReadObject();
		}

		// public void mouseWheelMoved(MouseWheelEvent e) {
		//		
		// int lockedX = (int)e.getPoint().getX();
		// if(e.getWheelRotation() > 0){
		// lockedX = (int)getScreenWidth() - lockedX;
		// }
		//		
		// float pointerBp = views.trackToBp(lockedX);
		// float pointerRelative = views.trackToRelative(lockedX);
		//		
		// long startBp = views.getBpRegion().start;
		// long endBp = views.getBpRegion().end;
		//		
		// float width = endBp - startBp;
		// width *= Math.pow(ZOOM_FACTOR, e.getWheelRotation());
		//		
		//		
		// startBp = (long)(pointerBp - width * pointerRelative);
		// endBp = (long)(pointerBp + width * (1 - pointerRelative));
		//		
		// if(startBp < 0){
		// endBp += -startBp;
		// startBp = 0;
		// }
		//		
		// //System.out.println(startBp + ", " + endBp);
		//		
		// views.setBpRegion(new Region(startBp, endBp));
		// }

		public void redraw() {
			this.datasetChanged(new DatasetChangeEvent(this, null));
		}

		public Collection<View> getViews() {
			return views;
		}
	}
