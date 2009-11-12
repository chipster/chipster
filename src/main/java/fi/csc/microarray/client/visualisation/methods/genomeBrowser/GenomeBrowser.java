package fi.csc.microarray.client.visualisation.methods.genomeBrowser;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
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

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.StraightforwardFileParser;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.TreeThread;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.BedInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.CasavaReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ConstantRefGeneInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.CytobandReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ElandReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.FastaFsfReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.FastaReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.EmptyTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.GeneTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.IntensityTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.SeparatorTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.SeqBlockTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.SeqTrack;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.track.GeneTrack.PartColor;


public class GenomeBrowser extends Plot 
implements ChartMouseListener, Cloneable, Serializable{ //, MouseWheelListener {	

	/** The cell information text is drawn with this font. */
	private Font descriptionFont;

	private List<View> views = new LinkedList<View>();

	private View overview;;

	public GenomeBrowser(){

		overview = new HorizontalView(this, false, false, true);
		//overview = new HorizontalView(this, false, false, true);
		
		File cytobandFile = new File("cytoband_hg17_sorted.txt");

		CytobandTrack overviewCytobands = new CytobandTrack(overview, cytobandFile, 
				StraightforwardFileParser.class, new CytobandReadInstructions(), false);

		overview.addTrack(overviewCytobands);
		overviewCytobands.initializeListener();

		//overview.addTrack(new RulerTrack(overview));

		overview.margin = 5;

		this.views.add(overview);

		View dataView = null;

		//Horizontal or circular?
		if(true){

			dataView = new HorizontalView(this, true, true, false);

		} else {

			dataView = new CircularView(this, true, true, false);
			dataView.margin = 20;
			dataView.addTrack(new EmptyTrack(dataView, 30));
		}

		CytobandTrack cytobands = new CytobandTrack(dataView, cytobandFile, 
				StraightforwardFileParser.class, new CytobandReadInstructions(), true);
		
		dataView.addTrack(cytobands);
		cytobands.initializeListener();
		
//		Track annotationOverview = new BlockTrack(dataView, new File("chr1_uscs.fsf"), 
//				TreeThread.class, new ConstantRefGeneInstructions(), PartColor.CDS.c.darker(), 
//				10000000, Long.MAX_VALUE);
		
		Track annotationOverview = new IntensityTrack(dataView, new File("chr1_uscs.fsf"), 
				TreeThread.class, new ConstantRefGeneInstructions(), PartColor.CDS.c.darker(), 
				10000000);
		
		dataView.addTrack(annotationOverview);
		annotationOverview.initializeListener();
		
		Track annotation = new GeneTrack(dataView, new File("chr1_uscs.fsf"), 
				TreeThread.class, new ConstantRefGeneInstructions(), Color.DARK_GRAY, 10000000);

		dataView.addTrack(annotation);
		annotation.initializeListener();
		
		dataView.addTrack(new SeparatorTrack(dataView));
		
// Strand aware conciser not implemented yet		
//		Track annotationOverviewReversed = new BlockTrack(dataView, new File("chr1_uscs.fsf"), 
//				TreeThread.class, new ConstantRefGeneInstructions(), PartColor.CDS.c.darker(), 
//				10000000, Long.MAX_VALUE);
//		
//		annotationOverviewReversed.setReverseStrand(true);
//		dataView.addTrack(annotationOverviewReversed);
//		annotationOverviewReversed.initializeListener();
		
		Track annotationReversed = new GeneTrack(dataView, new File("chr1_uscs.fsf"), 
				TreeThread.class, new ConstantRefGeneInstructions(), Color.DARK_GRAY, 10000000);

		annotationReversed.setReverseStrand(true);
		dataView.addTrack(annotationReversed);
		annotationReversed.initializeListener();
		
		//Eland export
		if(true){

			File userData = new File("eland_export.fsf");
			ReadInstructions<Float> userDataInstructions = new ElandReadInstructions();


			IntensityTrack reads2 = new IntensityTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.gray, 1000000);

			dataView.addTrack(reads2);
			reads2.initializeListener();

			SeqBlockTrack reads = new SeqBlockTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.RED, 0, 1000000);

			dataView.addTrack(reads);
			reads.initializeListener();
			
			dataView.addTrack(new SeparatorTrack(dataView));
			
			
			//Refernce sequence			
			File seqFile = new File("chr1.fsf");

			SeqTrack seq = new SeqTrack(dataView, seqFile, 
					TreeThread.class, new FastaFsfReadInstructions(), 800);

			dataView.addTrack(seq);
			seq.initializeListener();
			
			dataView.addTrack(new SeparatorTrack(dataView));
			
			SeqBlockTrack readsReversed = new SeqBlockTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.RED, 0, 1000000);

			readsReversed.setReverseStrand(true);
			dataView.addTrack(readsReversed);
			readsReversed.initializeListener();
			
		}

		
		//Casava file, no peaks
		if(false){

			File userData = new File("casava_chr1_filtered.fsf");
			ReadInstructions<Float> userDataInstructions = new CasavaReadInstructions();


			IntensityTrack reads2 = new IntensityTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.gray, 1000000);

			dataView.addTrack(reads2);
			reads2.initializeListener();

			SeqBlockTrack reads = new SeqBlockTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.RED, 0, 1000000);

			dataView.addTrack(reads);
			reads.initializeListener();
			
			dataView.addTrack(new SeparatorTrack(dataView));
			
			
			//Original sequence
			
			File seqFile = new File("chr1.fa");
			try {
				SeqTrack seq = new SeqTrack(dataView, seqFile, 
						TreeThread.class, new FastaReadInstructions(seqFile), 400);
				
				dataView.addTrack(seq);
				seq.initializeListener();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			dataView.addTrack(new SeparatorTrack(dataView));
			
			SeqBlockTrack readsReversed = new SeqBlockTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.RED, 0, 1000000);

			readsReversed.setReverseStrand(true);
			dataView.addTrack(readsReversed);
			readsReversed.initializeListener();
			
		}
		
		
		//BED file
		if(false){
			File userData = new File("bcell.fsf");
			ReadInstructions<Float> userDataInstructions = new BedInstructions();


			IntensityTrack reads2 = new IntensityTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.gray, 1000000);

			dataView.addTrack(reads2);
			reads2.initializeListener();

			SeqBlockTrack reads = new SeqBlockTrack(dataView, userData, 
					TreeThread.class, userDataInstructions, Color.RED, 0, 1000000);

			dataView.addTrack(reads);
			reads.initializeListener();
			
		}
		
//		IntensityTrack line = new IntensityTrack(dataView, userData, 
//				TreeThread.class, userDataInstructions, 1024*1024);
//
//		dataView.addTrack(line);
//		line.initializeListener();


		dataView.addTrack(new RulerTrack(dataView));
		
		this.views.add(dataView);    	    	    	        

		dataView.addRegionListener(new RegionListener(){
			public void RegionChanged(Region bpRegion) {
				overview.highlight = bpRegion;				
			}    		
		});
	}



	public String getPlotType() {
		return "GeneBrowser";
	}


	public Font getDescriptionFont() {
		return this.descriptionFont;
	}

	/**
	 * Draws the plot on a Java2D graphics device (such as the screen or 
	 * a printer).
	 *
	 * @param g2  the graphics device.
	 * @param area  the area within which the plot should be drawn.
	 * @param anchor  the anchor point (<code>null</code> permitted).
	 * @param parentState  the state from the parent plot, if there is one
	 * (<code>null</code> permitted.)
	 * @param info  collects info about the drawing (<code>null</code> permitted).
	 * @throws NullPointerException  if g2 or area is null.
	 */ 
	public void draw(java.awt.Graphics2D g2, java.awt.geom.Rectangle2D area, 
			java.awt.geom.Point2D anchor, PlotState parentState, 
			PlotRenderingInfo info) {

		// adjust for insets...
		//this.getInsets().trim(area);

		if (info != null) {
			info.setPlotArea(area);
			info.setDataArea(area);
		}

		//this.setBackgroundPaint(Color.black);

		drawBackground(g2, area);
		drawOutline(g2, area);

		Shape savedClip = g2.getClip();
		g2.clip(area);
		/*
        Composite originalComposite = g2.getComposite();
        g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 
                getForegroundAlpha()));
		 */
		drawPlot(g2, area, info);

		Rectangle viewArea = (Rectangle) area.getBounds().clone();

		//Horizontal or vertical split
		if(true){

			float[] viewHeights = new float[] { 0.2f, 0.8f };
			int origHeight = viewArea.height;

			for(int i = 0 ; i < views.size(); i++){
				if( i > 0){
					viewArea.y += viewArea.height;
				}
				viewArea.height = (int)(area.getBounds().getHeight() * viewHeights[i]);
				g2.setClip(viewArea);
				views.get(i).drawView(g2, false);
			}
		}else {
			float[] viewWidths = new float[] { 0.05f, 0.95f };
			Rectangle lastArea = null;

			for(int i = 0 ; i < views.size(); i++){
				if( lastArea != null){
					viewArea.x = lastArea.x + lastArea.width;
					viewArea.y = lastArea.y;
					viewArea.height = lastArea.height;
				}
				g2.setColor(Color.black);
				viewArea.width = (int)(area.getBounds().getWidth() * viewWidths[i]);
				lastArea = (Rectangle)(viewArea.clone());

				View view = views.get(i);       

				if(view instanceof VerticalView){
					viewArea.grow(0, -view.margin);
				} else if (view instanceof HorizontalView){
					viewArea.grow(-view.margin, 0);
				}

				g2.setClip(savedClip);
				g2.drawLine(viewArea.x - 1, 0, viewArea.x - 1, viewArea.height);

				g2.setClip(viewArea);
				view.drawView(g2, false);
			}
		}
		g2.setClip(savedClip);
		//   g2.setComposite(originalComposite);

		drawOutline(g2, area);
	}

	private int screenWidth;

	public ChartPanel chartPanel;

	private void drawPlot(Graphics2D g2, java.awt.geom.Rectangle2D area, 
			PlotRenderingInfo info) {

		//if(view.getBpRegion()  != null){
		//    		g2.drawString("Showing: " + 
		//    				Utils.toHumanReadable((view.getBpRegion().end - view.getBpRegion().start) * 51)  + 
		//    				", RAM footprint: " + Utils.toHumanReadable(Runtime.getRuntime().totalMemory()), 10, 20);
		//}
		//    	if(lastStatus != null){
		//    		//g2.fillRect(0, buf.getHeight() - 30, buf.getWidth() - (int)lastStatus.areaRequestCount, 10);
		//    		g2.fillRect(0, buf.getHeight() - 5, buf.getWidth() - (int)lastStatus.fileRequestCount, 5);
		//    		//g2.fillRect(0, buf.getHeight() - 10, buf.getWidth() - (int)lastStatus.fileResultCount, 10);
		//    	}
	}

	private float getScreenWidth() {
		if(screenWidth != 0){
			return screenWidth;
		} else {
			return 600;
		}
	}


	/**
	 * Implements the ChartMouseListener interface. This 
	 * method does nothing.
	 *
	 * @param event  the mouse event.
	 */ 
	public void chartMouseMoved(ChartMouseEvent event) {
	}

	public void chartMouseClicked(ChartMouseEvent e) {
	}

	/**
	 * Tests this plot for equality with an arbitrary object.  Note that the 
	 * plot's dataset is NOT included in the test for equality.
	 *
	 * @param obj  the object to test against (<code>null</code> permitted).
	 *
	 * @return <code>true</code> or <code>false</code>.
	 */
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		if (!(obj instanceof GenomeBrowser)) {
			return false;
		}
		GenomeBrowser that = (GenomeBrowser) obj;

		//        if (!ObjectUtilities.equal(this.urlGenerator, that.urlGenerator)) {
		//            return false;
		//        }

		if (!ObjectUtilities.equal(this.descriptionFont, 
				that.descriptionFont)) {
			return false;
		}

		// can't find any difference...
		return true;
	}

	/**
	 * Provides serialization support.
	 *
	 * @param stream  the output stream.
	 *
	 * @throws IOException  if there is an I/O error.
	 * @throws NullPointerException  if stream is null.
	 */
	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.defaultWriteObject();
	}

	/**
	 * Provides serialization support.
	 *
	 * @param stream  the input stream.
	 *
	 * @throws IOException  if there is an I/O error.
	 * @throws ClassNotFoundException  if there is a classpath problem.
	 * @throws NullPointerException  if stream is null.
	 */
	private void readObject(ObjectInputStream stream) 
	throws IOException, ClassNotFoundException {
		stream.defaultReadObject();
	}





	//	public void mouseWheelMoved(MouseWheelEvent e) {	
	//		
	//		int lockedX = (int)e.getPoint().getX();
	//		if(e.getWheelRotation() > 0){
	//			lockedX = (int)getScreenWidth() - lockedX;
	//		}
	//		
	//		float pointerBp = views.trackToBp(lockedX);
	//		float pointerRelative = views.trackToRelative(lockedX);
	//		
	//		long startBp = views.getBpRegion().start;
	//		long endBp = views.getBpRegion().end;
	//		
	//		float width = endBp - startBp;
	//		width *= Math.pow(ZOOM_FACTOR, e.getWheelRotation());			
	//		
	//		
	//		startBp = (long)(pointerBp - width * pointerRelative);
	//		endBp = (long)(pointerBp + width * (1 - pointerRelative));
	//		
	//		if(startBp < 0){
	//			endBp += -startBp;
	//			startBp = 0;
	//		}
	//		
	//		//System.out.println(startBp + ", " + endBp);
	//		
	//		views.setBpRegion(new Region(startBp, endBp));
	//	}

	public void redraw(){

		this.datasetChanged(new DatasetChangeEvent(this, null));
	}


	public Collection<View> getViews() {
		return views;
	}
}

