package fi.csc.microarray.client.visualisation.methods;

import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRendererState;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleEdge;
import org.jfree.util.ShapeUtilities;
import org.jfree.util.UnitType;

import fi.csc.microarray.client.visualisation.methods.Scatterplot.DataItem2D;

public class PositionRecordingRenderer extends StandardXYItemRenderer {

	protected List<DataItem2D> allItems = new LinkedList<DataItem2D>();
	protected Set<DataItem2D> selectedItems = new HashSet<DataItem2D>();
	protected ChartPanel chartPanel;

	public PositionRecordingRenderer(int shapes, List<DataItem2D> allItems, Set<DataItem2D> selectedItems, ChartPanel chartPanel) {
		super(shapes);
		this.allItems = allItems;
		this.selectedItems = selectedItems;
		this.chartPanel = chartPanel;
	}

	@Override
	public void drawItem(Graphics2D g2, XYItemRendererState state, Rectangle2D dataArea, PlotRenderingInfo info, XYPlot plot, ValueAxis domainAxis, ValueAxis rangeAxis, XYDataset dataset, int series, int item, CrosshairState crosshairState, int pass) {

		// this method body has been ripped from JFreeChart, and modified to
		// work in a derived class. it does not support every feature the
		// original did

		boolean itemVisible = getItemVisible(series, item);

		// setup for collecting optional entity info...
		Shape entityArea = null;
		EntityCollection entities = null;
		if (info != null) {
			entities = info.getOwner().getEntityCollection();
		}

		PlotOrientation orientation = plot.getOrientation();
		Paint paint = getItemPaint(series, item);
		Stroke seriesStroke = getItemStroke(series, item);
		g2.setPaint(paint);
		g2.setStroke(seriesStroke);

		// get the data point...
		double x1 = dataset.getXValue(series, item);
		double y1 = dataset.getYValue(series, item);
		if (Double.isNaN(x1) || Double.isNaN(y1)) {
			itemVisible = false;
		}

		RectangleEdge xAxisLocation = plot.getDomainAxisEdge();
		RectangleEdge yAxisLocation = plot.getRangeAxisEdge();
		double transX1 = domainAxis.valueToJava2D(x1, dataArea, xAxisLocation);
		double transY1 = rangeAxis.valueToJava2D(y1, dataArea, yAxisLocation);

		// logger.debug("x1: " + x1 + "y1: " + y1 + "transX1: " + transX1 +
		// "transY1: " + transY1);

		if (getPlotLines()) {
			if (true/* this.drawSeriesLineAsPath */) {
				State s = (State) state;
				if (true/* s.getSeriesIndex() != series */) {
					// we are starting a new series path
					/*
					 * s.seriesPath().reset(); s.lastPointGood = false;
					 * s.setSeriesIndex(series);
					 */
				}

				// update path to reflect latest point
				if (itemVisible && !Double.isNaN(transX1) && !Double.isNaN(transY1)) {
					float x = (float) transX1;
					float y = (float) transY1;
					if (orientation == PlotOrientation.HORIZONTAL) {
						x = (float) transY1;
						y = (float) transX1;
					}
					if (s.isLastPointGood()) {
						// TODO: check threshold
						s.seriesPath.lineTo(x, y);
					} else {
						s.seriesPath.moveTo(x, y);
					}
					s.setLastPointGood(true);
				} else {
					s.setLastPointGood(false);
				}
				if (item == dataset.getItemCount(series) - 1) {
					if (true/* s.seriesIndex == series */) {
						// draw path
						g2.setStroke(getSeriesStroke(series));
						g2.setPaint(getSeriesPaint(series));
						g2.draw(s.seriesPath);
					}
				}
			}

			else if (item != 0 && itemVisible) {
				// get the previous data point...
				double x0 = dataset.getXValue(series, item - 1);
				double y0 = dataset.getYValue(series, item - 1);
				if (!Double.isNaN(x0) && !Double.isNaN(y0)) {
					boolean drawLine = true;
					if (getPlotDiscontinuous()) {
						// only draw a line if the gap between the current and
						// previous data point is within the threshold
						int numX = dataset.getItemCount(series);
						double minX = dataset.getXValue(series, 0);
						double maxX = dataset.getXValue(series, numX - 1);
						if (this.getGapThresholdType() == UnitType.ABSOLUTE) {
							drawLine = Math.abs(x1 - x0) <= this.getGapThreshold();
						} else {
							drawLine = Math.abs(x1 - x0) <= ((maxX - minX) / numX * getGapThreshold());
						}
					}
					if (drawLine) {
						double transX0 = domainAxis.valueToJava2D(x0, dataArea, xAxisLocation);
						double transY0 = rangeAxis.valueToJava2D(y0, dataArea, yAxisLocation);

						// only draw if we have good values
						if (Double.isNaN(transX0) || Double.isNaN(transY0) || Double.isNaN(transX1) || Double.isNaN(transY1)) {
							return;
						}

						if (orientation == PlotOrientation.HORIZONTAL) {
							state.workingLine.setLine(transY0, transX0, transY1, transX1);
						} else if (orientation == PlotOrientation.VERTICAL) {
							state.workingLine.setLine(transX0, transY0, transX1, transY1);
						}

						if (state.workingLine.intersects(dataArea)) {
							g2.draw(state.workingLine);
						}
					}
				}
			}
		}

		// we needed to get this far even for invisible items, to ensure that
		// seriesPath updates happened, but now there is nothing more we need
		// to do for non-visible items...
		if (!itemVisible) {
			return;
		}

		if (getBaseShapesVisible()) {

			Shape shape = getItemShape(series, item);
			if (orientation == PlotOrientation.HORIZONTAL) {
				shape = ShapeUtilities.createTranslatedShape(shape, transY1, transX1);
			} else if (orientation == PlotOrientation.VERTICAL) {
				shape = ShapeUtilities.createTranslatedShape(shape, transX1, transY1);
			}
			if (shape.intersects(dataArea)) {
				if (getItemShapeFilled(series, item)) {
					g2.fill(shape);
				} else {
					g2.draw(shape);
				}

				boolean selected = false;
				for (DataItem2D selectedItem : selectedItems) {
//					System.out.println("drawing " + item + " " + series);
//					System.out.println("comparing to currently selected " + selectedItem.getName() + " " + selectedItem.getIndex() + " " + selectedItem.getSeries());
					if (selectedItem.getIndex() == item && selectedItem.getSeries() == series) {
						selected = true;
						break;
					}
				}

				if (selected) {
					Rectangle rect = shape.getBounds();
					rect.grow(2, 2);
					g2.draw(rect);
				}

				// store names and coordinates for point selection
				Rectangle rect = new Rectangle();
				rect.setSize(shape.getBounds().getSize());// TODO This should
				// be scaled too
				rect.setLocation(chartPanel.translateJava2DToScreen(shape.getBounds().getLocation()));

				allItems.get(item).setBounds(rect);

				// logger.debug(chartPanel.translateJava2DToScreen(shape.getBounds().getLocation()));
				// logger.debug("drew item " + itemNames.get(item) + " to " +
				// shape.getBounds().x + "," + shape.getBounds().y + ", with h="
				// + shape.getBounds().height + ", w=" +
				// shape.getBounds().width);
			}
			entityArea = shape;

		}

		if (getPlotImages()) {
			Image image = getImage(plot, series, item, transX1, transY1);
			if (image != null) {
				Point hotspot = getImageHotspot(plot, series, item, transX1, transY1, image);
				g2.drawImage(image, (int) (transX1 - hotspot.getX()), (int) (transY1 - hotspot.getY()), null);
				entityArea = new Rectangle2D.Double(transX1 - hotspot.getX(), transY1 - hotspot.getY(), image.getWidth(null), image.getHeight(null));
			}

		}

		// draw the item label if there is one...
		if (isItemLabelVisible(series, item)) {
			double xx = transX1;
			double yy = transY1;
			if (orientation == PlotOrientation.HORIZONTAL) {
				xx = transY1;
				yy = transX1;
			}
			drawItemLabel(g2, orientation, dataset, series, item, xx, yy, (y1 < 0.0));
		}

		/*
		 * int domainAxisIndex = plot.getDomainAxisIndex(domainAxis); int
		 * rangeAxisIndex = plot.getRangeAxisIndex(rangeAxis);
		 * updateCrosshairValues(crosshairState, x1, y1, domainAxisIndex,
		 * rangeAxisIndex, transX1, transY1, orientation);
		 */

		// add an entity for the item...
		if (entities != null) {
			addEntity(entities, entityArea, dataset, series, item, transX1, transY1);
		}
	}

}
