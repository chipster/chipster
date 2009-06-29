package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.HCPlot;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.OverlayLayout;
import org.jfree.ui.RectangleEdge;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.Visualisation;

/**
 * JFreeChart handles area selections with zoom whereas we wan't selection to hapen. This
 * class wraps and extends ChartPanel functionality. 
 * 
 * @author klemela
 *
 */
public class SelectableChartPanel extends JPanel implements MouseListener, MouseMotionListener {

	private TransparentPanel transparentPanel;
	private ChartPanel chartPanel;
	
	private SelectionChangeListener selectionListener;
	private Rectangle2D.Double selection;
	
	//Maybe in future selectable charts have also the zoom functionality
	private final boolean isZoomable = false;
	
	//This class has it's own overlay panel for drawing selection rectangles to
	//make horizontal selection possible in category plots. With
	//XYPlot it's also possible to use JFC:s zoom outline
	private final boolean useZoomOutline = false;
	
	public interface SelectionChangeListener{
		/**
		 * Data points inside given rectangle should be added to selection if they aren't
		 * there already and removed if they are. If the rectangle is null, selection should be 
		 * cleared. This makes it possible to handle all situations with single logic. If selection
		 * is made without pressing control key, this method is called first with null argument to 
		 * clear the selection and again with new selection rectangle.
		 * 
		 * @param newSelection
		 */
		public void selectionChanged(Rectangle2D.Double newSelection);
	}
	
	public SelectableChartPanel(JFreeChart chart, SelectionChangeListener selectionListener) {
		this(chart, selectionListener, true);
	}
	
	public SelectableChartPanel(JFreeChart chart, SelectionChangeListener selectionListener, boolean scalable) {
		super(new OverlayLayout());
		
		this.selectionListener = selectionListener;
		
		if(scalable){
			chartPanel = Visualisation.makePanel(chart);
		} else {
			chartPanel = Visualisation.makenNonScalablePanel(chart);
		}
		transparentPanel = new TransparentPanel();
		
		this.add(transparentPanel);
		this.add(chartPanel);
		
		chartPanel.addMouseListener(this);
		chartPanel.addMouseMotionListener(this);
		chartPanel.setMouseZoomable(isZoomable);
	}
	
	public Rectangle.Double translateToChart(Rectangle mouseCoords){
		
		//Screen coordinates are integers
		Point corner1 = new Point((int)mouseCoords.getMinX(), (int)mouseCoords.getMinY());
		Point corner2 = new Point((int)mouseCoords.getMaxX(), (int)mouseCoords.getMaxY());
		
		Point.Double translated1 = translateToChart(corner1);
		Point.Double translated2 = translateToChart(corner2);
		
		return createOrderedRectangle(translated1, translated2);				
	}


	public Point2D.Double translateToChart(Point mouseCoords)
	{
		Insets insets = getInsets();
		int mouseX = (int) ((mouseCoords.getX() - insets.left) / chartPanel.getScaleX());
		int mouseY = (int) ((mouseCoords.getY() - insets.top) / chartPanel.getScaleY());

		Point2D p = chartPanel.translateScreenToJava2D(new Point(mouseX, mouseY));	
		
		Plot plot = chartPanel.getChart().getPlot();
		
		if (plot instanceof XYPlot) { 
						
			XYPlot xyPlot = (XYPlot) plot;
			
			ChartRenderingInfo info = chartPanel.getChartRenderingInfo();
			Rectangle2D dataArea = info.getPlotInfo().getDataArea();

			ValueAxis domainAxis = xyPlot.getDomainAxis();
			RectangleEdge domainAxisEdge = xyPlot.getDomainAxisEdge();
			ValueAxis rangeAxis = xyPlot.getRangeAxis();
			RectangleEdge rangeAxisEdge = xyPlot.getRangeAxisEdge();						
			
			/*multiplication by getScale to take care of the scaling of Java2D Graphics from Java2D
			 * coordinates to screen coordinates. Function java2DToValue is used to take care of the
			 * internal scaling of JFreeChart which scales data values to Java2D coordinates 
			*/
			
			double chartX = domainAxis.java2DToValue(p.getX() *	chartPanel.getScaleX(), 
					dataArea, domainAxisEdge);
			
			double chartY = rangeAxis.java2DToValue(p.getY() * chartPanel.getScaleY(), 
					dataArea, rangeAxisEdge);
			
			return new Point2D.Double(chartX, chartY);
			
		} else if (plot instanceof CategoryPlot) {		
			
			CategoryPlot catPlot = (CategoryPlot) plot;
			
			ChartRenderingInfo info = chartPanel.getChartRenderingInfo();
			Rectangle2D dataArea = info.getPlotInfo().getDataArea();

			CategoryAxis domainAxis = catPlot.getDomainAxis();
			RectangleEdge domainAxisEdge = catPlot.getDomainAxisEdge();
			int categoryCount = catPlot.getCategories().size();
			ValueAxis rangeAxis = catPlot.getRangeAxis();
			RectangleEdge rangeAxisEdge = catPlot.getRangeAxisEdge();

			//See comment from getSelectionRectangle method for x coordinate translation
			
			double firstCategoryX = domainAxis.getCategoryJava2DCoordinate(
					catPlot.getDomainGridlinePosition(), 0, categoryCount, dataArea, 
					domainAxisEdge);
			
			double lastCategoryX  = domainAxis.getCategoryJava2DCoordinate(
					catPlot.getDomainGridlinePosition(), categoryCount - 1, categoryCount, 
					dataArea, domainAxisEdge);
			
			/*multiplication by getScale to take care of the scaling of Java2D Graphics from Java2D
			 * coordinates to screen coordinates. Function java2DToValue is used to take care of the
			 * internal scaling of JFreeChart which scales data values to Java2D coordinates 
			*/
			
			double scaledX = p.getX() * chartPanel.getScaleX();
			double chartY = rangeAxis.java2DToValue(p.getY() * chartPanel.getScaleY(),
					dataArea, rangeAxisEdge);
			
			//Category axis doesn't have integer values, so we calculate it 
			double relativeX = (scaledX - firstCategoryX) / (lastCategoryX - firstCategoryX) * 
				(categoryCount - 1);
									
						
			return new Point2D.Double(relativeX, chartY);
		} else if (plot instanceof HCPlot) {		
			//HCPlot hcPlot = (HCPlot)plot;
			
	        Insets chartInsets = getInsets();
	        int x = (int) ((p.getX() - chartInsets.left) * chartPanel.getScaleX());
	        int y = (int) ((p.getY() - chartInsets.top) * chartPanel.getScaleY());
											
			return new Point2D.Double(x, y);
		} else {
			throw new UnsupportedOperationException("Only Category, XY and HC plots are supported until now");			
		}
	}
	
	/**
	 * Gives two corners of selection rectangle relative to original data values. If the 
	 * chart is categoryChart, the x values are scaled to category indexes so that click on first
	 * category returns 0, click on last return (categoryCount - 1) and other values are interpolated
	 * or extrapolated from these. If the selection was made with a single click, the width and 
	 * height have value of 0;
	 * 
	 * @return Selection rectangle, null if selection should be cleared
	 */
	public Rectangle2D.Double getSelectionRectangle(){
		return selection;
	}
	
	public void mouseClicked(MouseEvent e) {
		if (!e.isControlDown()) {
			setSelection(null);
		}
		
		//Rectangle made by single click is grown by couple pixels to make clicking easier.
		//Growing should be done by screen pixels, and this is the latest point for that for now.
		//If accurate single click detection is needed some day later, this should be moved to 
		//SelectionChangeListeners and implement needed methods or parameters.
		
		Rectangle rect = new Rectangle((int)e.getPoint().getX(), (int)e.getPoint().getY(), 0 ,0);
		rect.grow(3, 3);
		
		Rectangle.Double translatedRect = translateToChart(rect);
				
		setSelection(translatedRect);
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	private Point startCoords;
	private boolean isDragged = false;

	
	public void mousePressed(MouseEvent e)
	{
		if (useZoomOutline){

			chartPanel.setMouseZoomable(false,true);
			chartPanel.mousePressed(e);			
			chartPanel.setMouseZoomable(true,false);

		}
		if(isZoomable){
			chartPanel.mousePressed(e);
		}

		startCoords = e.getPoint();
		isDragged = false;
	}
	
	public void mouseReleased(MouseEvent e) {
		
		if (useZoomOutline){
			chartPanel.setMouseZoomable(false,true);
			chartPanel.mouseReleased(e);
			chartPanel.setMouseZoomable(true,false);

		}
		if(isZoomable){
			chartPanel.mouseReleased(e);
		}
		
		// Hide the selection rectangle
		transparentPanel.setArea(null);
		transparentPanel.repaint();

		if (isDragged) {
			if (!e.isControlDown()) {
				setSelection(null);
			}
			
			Point2D.Double translatedStart = translateToChart(startCoords);
			Point2D.Double translatedEnd  = translateToChart(e.getPoint());
					
			setSelection(createOrderedRectangle(translatedStart, translatedEnd));
		}
	}

	public void mouseDragged(MouseEvent e) {
		isDragged = true;			
		transparentPanel.setArea(createOrderedRectangle(startCoords, e.getPoint()).getBounds());
		transparentPanel.repaint();
	}

	public void mouseMoved(MouseEvent e) {
	}
	
	/**
	 * Creates rectangle, in which the first coordinate pair is upper left corner
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	private Rectangle2D.Double createOrderedRectangle(Point2D p1, Point2D p2){
		double x = p1.getX() < p2.getX() ? p1.getX() : p2.getX();
		double y = p1.getY() < p2.getY() ? p1.getY() : p2.getY();
		double w = Math.abs(p2.getX() - p1.getX());
		double h = Math.abs(p2.getY() - p1.getY());
		
		return new Rectangle2D.Double(x, y, w, h);
	}
	
	private void setSelection(Rectangle2D.Double selection){
		this.selection = selection;
		if(selectionListener != null){						
			
			selectionListener.selectionChanged(this.selection);
		}
	}
	
	protected class TransparentPanel extends JPanel {

		private Rectangle2D area;

		public TransparentPanel() {
			this.setOpaque(false);
		}

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponents(g);

			if (area != null) {
				Graphics2D g2d = (Graphics2D) g;
				g2d.setColor(Color.DARK_GRAY);
				g2d.setStroke(VisualConstants.dashLine);

				g2d.draw(area);
			}
		}

		protected void setArea(Rectangle area) {
			this.area = area;
		}
	}

	public ChartPanel getChartPanel() {
		return chartPanel;
	}
}
