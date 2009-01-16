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
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.OverlayLayout;
import org.jfree.ui.RectangleEdge;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.Visualisation;

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
		public void selectionChanged(Rectangle2D.Double newSelection);
	}
	
	public SelectableChartPanel(JFreeChart chart, SelectionChangeListener selectionListener) {
		super(new OverlayLayout());
		
		this.selectionListener = selectionListener;
		
		chartPanel = Visualisation.makePanel(chart);		
		transparentPanel = new TransparentPanel();
		
		this.add(transparentPanel);
		this.add(chartPanel);
		
		chartPanel.addMouseListener(this);
		chartPanel.addMouseMotionListener(this);
		chartPanel.setMouseZoomable(isZoomable);
	}


	public Point2D.Double translateToChart(Point mouseCoords)
	{
		Insets insets = getInsets();
		int mouseX = (int) ((mouseCoords.getX() - insets.left) / chartPanel.getScaleX());
		int mouseY = (int) ((mouseCoords.getY() - insets.top) / chartPanel.getScaleY());

		Point2D p = chartPanel.translateScreenToJava2D(new Point(mouseX, mouseY));	
		
		Plot plot = chartPanel.getChart().getPlot();
		
		if (plot instanceof XYPlot) {
			//TODO Replace hackish PositionRecordingRenderer from Scatterplot and Volcano plot 
			// with this one. These coordinates are in original scale, and can be compared directly 
			// to original data, so collection of rendered coordinates isn't needed anymore 
						
			XYPlot xyPlot = (XYPlot) plot;
			
			ChartRenderingInfo info = chartPanel.getChartRenderingInfo();
			Rectangle2D dataArea = info.getPlotInfo().getDataArea();

			ValueAxis domainAxis = xyPlot.getDomainAxis();
			RectangleEdge domainAxisEdge = xyPlot.getDomainAxisEdge();
			ValueAxis rangeAxis = xyPlot.getRangeAxis();
			RectangleEdge rangeAxisEdge = xyPlot.getRangeAxisEdge();

			
			double chartX = domainAxis.java2DToValue(p.getX(), dataArea, domainAxisEdge);
			double chartY = rangeAxis.java2DToValue(p.getY(), dataArea, rangeAxisEdge);
			
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
			
			
			double relativeX = (p.getX() - firstCategoryX) / (lastCategoryX - firstCategoryX) * 
				(categoryCount - 1);
						
			
			double chartY = rangeAxis.java2DToValue(p.getY(), dataArea, rangeAxisEdge);
						
			return new Point2D.Double(relativeX, chartY);
		} else {
			throw new UnsupportedOperationException("Only Category and XY plots are supported");			
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
		
		Point2D p = translateToChart(e.getPoint());		
		setSelection(new Rectangle2D.Double(p.getX(), p.getY(), 0, 0));
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
			System.out.println("mouse pressed");
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

		//TODO this should be done in the selection listeners
		// Draw the selection frames for the dataItems
		//chartPanel.getChart().fireChartChanged();
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
}
