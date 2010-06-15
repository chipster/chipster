package fi.csc.microarray.client.dataview;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.JViewport;
import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;
import org.jgraph.event.GraphSelectionEvent;
import org.jgraph.event.GraphSelectionListener;
import org.jgraph.graph.BasicMarqueeHandler;
import org.jgraph.graph.CellViewFactory;
import org.jgraph.graph.DefaultCellViewFactory;
import org.jgraph.graph.DefaultGraphModel;
import org.jgraph.graph.GraphLayoutCache;
import org.jgraph.graph.GraphModel;
import org.jgraph.graph.VertexView;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.ToolBarComponentFactory;
import fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GraphRenderer;
import fi.csc.microarray.client.dataviews.vertexes.GraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GroupVertex;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.LinksChangedEvent;
import fi.csc.microarray.util.SwingTools;

/**
 * A GUI component for viewing the data in a graph presentation, where each dataset is represented by a graphical vertex (an instance of
 * GraphVertex). The graph is created by using the JGraph application framework.
 * 
 * @author Janne KÃ¤ki, Aleksi Kallio, Petri KlemelÃ¤
 * 
 */
public class GraphPanel extends JPanel implements ActionListener, PropertyChangeListener, DataChangeListener, AnimatorScrollable {

	public final float ZOOM_FACTOR = 1.2f;
	public final float ZOOM_IN_LIMIT = 1.0f;
	public final float ZOOM_OUT_LIMIT = 0.2f;

	private MicroarrayGraph graph = null;

	private ClientApplication application = Session.getSession().getApplication();
	private GraphModel model = new DefaultGraphModel(); // FIXME memory leak

	private JScrollPane graphScroller = null;

	private JToolBar buttonToolBar = null;
	private JButton zoomInButton;
	private JButton zoomOutButton;
	private JCheckBox autoZoomChecBbox;

	private JButton historyButton;

	private boolean internalSelection = false;

	private static final Logger logger = Logger.getLogger(GraphPanel.class);

	public void setInternalSelection(boolean internalSelection) {
		this.internalSelection = internalSelection;
	}

	private static final double DEFAULT_GRID_SIZE = 10.0;
	private static final double BIGGER_GRID_SIZE = 20.0;
	private static final double HUGE_GRID_SIZE = 40.0;

	/**
	 * Creates a new GraphPanel with default contents and appearance.
	 */
	public GraphPanel() {
		
		getGraph().getSelectionModel().addGraphSelectionListener(new WorkflowSelectionListener());
		getGraph().setMarqueeHandler(new BasicMarqueeHandler());

		// getGraph().setDebugGraphicsOptions(DebugGraphics.LOG_OPTION);
		// getGraph().setDoubleBuffered(false);

		this.setMinimumSize(new Dimension(0, 0));

		this.setLayout(new GridBagLayout());

		buttonToolBar = this.getButtonToolBar();
		graphScroller = this.getGraphScroller();

		this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.GRAPH_PANEL_HEIGHT));
		graphScroller.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.GRAPH_PANEL_HEIGHT));

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.SOUTH;
		c.gridy = 1;
		this.add(buttonToolBar, c);
		c.gridx = 0;
		c.gridy = 2;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.NORTH;
		c.weightx = 1.0;
		c.weighty = 1.0;
		this.add(graphScroller, c);

		// start listening
		application.addPropertyChangeListener(this);
		application.getDataManager().addDataChangeListener(this);

		graph.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));

		graph.setScale(graph.getScale() / ZOOM_FACTOR);

		// adds vertex renderer which adds the small '+' and '-' buttons to group
		VertexView.renderer = new GraphRenderer();
	}

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == zoomInButton) {
			autoZoomChecBbox.setSelected(false);
			this.zoomInToolUsed();

		} else if (source == zoomOutButton) {
			autoZoomChecBbox.setSelected(false);
			this.zoomOutToolUsed();

		} else if (source == autoZoomChecBbox) {
			autoZoom();
		} else if (source == historyButton) {
			application.showHistoryScreenFor(application.getSelectionManager().getSelectedDataBean());
		}
	}

	public void autoZoom() {
		if (autoZoomChecBbox.isSelected()) {
			Dimension dim = graph.getGraphSize();
			double xScale = graphScroller.getSize().getWidth() / (double) dim.width;
			double yScale = graphScroller.getSize().getHeight() / (double) dim.height;

			if (xScale < yScale) {
				this.setGraphScale(xScale);
			} else {
				this.setGraphScale(yScale);
			}

			this.updateGridSize();

			logger.debug("Graph size: " + dim.width + " , " + dim.height + " Scroller size: " + graphScroller.getSize().getWidth() + " , " + graphScroller.getSize().getHeight() + " Scales: " + xScale + " , " + yScale);

			graph.repaint();
		}
	}

	/**
	 * Will return (and if needed, create) the JGraph component that is the very heart of this GraphPanel. Initializes an empty graph and
	 * enforces differend kind of UI constraints on it (not editable, not connectable or disconnectable, not sizeable, etc., by user). Also
	 * makes the graph antialiased, which has a creat impact on the looks of this application :)
	 * 
	 * @return The graph component of this panel.
	 */
	public MicroarrayGraph getGraph() {
		if (graph == null) {

			// Sets missing constructor parameter for GraphLayout
			CellViewFactory factory = new DefaultCellViewFactory();
			boolean partial = true; // This must be true to allow group collapsing

			GraphLayoutCache cache = new GraphLayoutCache(model, factory, partial);

			graph = new MicroarrayGraph(model, cache, this);

			// Adds mouse listener which opens popup menu
			this.graph.addMouseListener(new MouseAdapter() {

				@Override
				public void mousePressed(MouseEvent e) {
					maybeShowPopup(e);
				}

				@Override
				public void mouseReleased(MouseEvent e) {
					maybeShowPopup(e);
				}

				public void mouseClicked(MouseEvent e) {

					logger.debug("mouseClicked");

					// Double click
					if (SwingUtilities.isLeftMouseButton(e) && e.getClickCount() > 1) {
						// Change the cursor back from marquee cross
						graph.setCursor(Cursor.getDefaultCursor());
						mouseButtonDoubleClicked(e);
					}
				}

				private void maybeShowPopup(MouseEvent e) {
					if (e.isPopupTrigger()) {
						showPopupMenu(e);
					}
				}
			});

			graph.setAntiAliased(true);
			graph.setBendable(false);
			graph.setConnectable(false);
			graph.setDisconnectable(false);
			graph.setEditable(false);
			graph.setMoveable(true);
			graph.setPortsVisible(false);
			graph.setSizeable(false);

			graph.setGridColor(Color.LIGHT_GRAY);
			graph.setGridSize(DEFAULT_GRID_SIZE);
			graph.setGridVisible(true);
			graph.setGridEnabled(true);

		}
		return graph;
	}

	private void mouseButtonDoubleClicked(MouseEvent e) {
		// Get the clicked cell
		Object cell = graph.getFirstCellForLocation(e.getX(), e.getY());

		logger.debug("Selected cell: " + cell);

		// Do not visualise collapsed group
		if (cell instanceof GraphVertex) {

			application.visualiseWithBestMethod(FrameType.MAIN);
		}
	}

	protected GraphModel getGraphModel() {
		return model;
	}

	/**
	 * Creates the surrounding scroller component for the graph. When called again, will simply return the existing scroller. The JGraph
	 * component will also be initialized from this method by calling getGraph(), if not previously created.
	 * 
	 * @return A JScrollPane which contains the graph component with appropriate margins.
	 */
	private JScrollPane getGraphScroller() {
		if (graphScroller == null) {
			graphScroller = new JScrollPane(this.getGraph());
			graphScroller.setBorder(BorderFactory.createEmptyBorder());
			graphScroller.setMinimumSize(new Dimension(0, 0));

			ScrollListener scrollListener = new ScrollListener();
			graphScroller.getHorizontalScrollBar().addAdjustmentListener(scrollListener);
			graphScroller.getVerticalScrollBar().addAdjustmentListener(scrollListener);
		}
		return graphScroller;
	}

	/**
	 * Repaints the graph after scroll to ensure the visibility of all vertexes
	 */
	private class ScrollListener implements AdjustmentListener {		
		// FIXME check if needed still with the new version of JGraph
		public void adjustmentValueChanged(AdjustmentEvent e) {
			GraphPanel.this.getGraph().repaint();
		}
	}

	/**
	 * When called for first time, creates button panel for the buttons of workflow -view.
	 * 
	 * @return The history panel of this GraphPanel.
	 */
	public JToolBar getButtonToolBar() {
		if (buttonToolBar == null) {

			zoomInButton = ToolBarComponentFactory.createButton(false, false);
			zoomInButton.setToolTipText("Zoom in");
			this.initialiseToolBarButton(zoomInButton);
			zoomInButton.setIcon(VisualConstants.ZOOM_IN_ICON);

			zoomOutButton = ToolBarComponentFactory.createButton(false, false);
			zoomOutButton.setToolTipText("Zoom out");
			this.initialiseToolBarButton(zoomOutButton);
			zoomOutButton.setIcon(VisualConstants.ZOOM_OUT_ICON);

			autoZoomChecBbox = ToolBarComponentFactory.createCheckBox("Fit");
			autoZoomChecBbox.setToolTipText("Scale workflow to show all datasets");
			autoZoomChecBbox.setSelected(true);
			this.initialiseToolBarButton(autoZoomChecBbox);

			historyButton = ToolBarComponentFactory.createButton(true, false);
			historyButton.setToolTipText("Show analyse history for the selected dataset");
			this.initialiseToolBarButton(historyButton);
			historyButton.setIcon(VisualConstants.GENERATE_HISTORY_ICON);
			historyButton.setEnabled(false);

			buttonToolBar = new JToolBar();
			buttonToolBar.setFloatable(false);
			buttonToolBar.setMinimumSize(new Dimension(0, 0));
			buttonToolBar.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);

			buttonToolBar.add(zoomInButton);
			buttonToolBar.add(zoomOutButton);
			buttonToolBar.add(autoZoomChecBbox);
			buttonToolBar.add(Box.createHorizontalGlue());
			buttonToolBar.add(historyButton);
		}
		return buttonToolBar;
	}

	private void initialiseToolBarButton(AbstractButton button) {
		button.addActionListener(this);
		button.setBorder(BorderFactory.createEmptyBorder(3, 3, 3, 3));
	}

	/**
	 * Centres the workflow view to given point.
	 * 
	 * @param p
	 *            Point must be given in the pixel coordinates of the graph. If this is done after scale of the graph, make sure that the
	 *            scaling has been done before calling this.
	 */
	public void pointToCenter(Point2D target) {
		Point newViewPos = new Point();

		newViewPos.x = (int) (target.getX() - graphScroller.getViewport().getWidth() / 2.0);
		newViewPos.y = (int) (target.getY() - graphScroller.getViewport().getHeight() / 2.0);

		this.setViewPosition(newViewPos);
	}

	/**
	 * Sets viewport position but checks that graph isn't moved outside the scrolling area
	 * 
	 * @param p
	 *            New View position
	 */
	public void setViewPosition(Point p) {
		Point2D p2 = graph.toScreen(p);

		JScrollPane scroller = this.getScroller();

		double newX = p2.getX();
		double newY = p2.getY();

		// magic numbers make some extra margin to avoid clipping
		int graphW = graph.getWidth();
		int graphH = graph.getHeight();
		int viewW = scroller.getViewport().getWidth();
		int viewH = scroller.getViewport().getHeight();

		// set the limits to prevent scrolling out of the area from the lower right corner
		int xLimit = graphW - viewW;
		int yLimit = graphH - viewH;

		// these are checked first, so that they are overriden if the graph is too small
		if (newX > xLimit) {
			newX = xLimit;
		}
		if (newY > yLimit) {
			newY = yLimit;
		}

		// set the limits to prevent scrolling out of the area from the upper left corner
		if (newX < 0) {
			newX = 0;
		}
		if (newY < 0) {
			newY = 0;
		}

		// apply chanches
		scroller.getViewport().setViewPosition(new Point((int) newX, (int) newY));
		graph.repaint();
	}

	public JScrollPane getScroller() {
		return this.graphScroller;
	}

	/**
	 * @param scale
	 *            This is set for the scale of graph if it is between zoom limits. Returns true if the scale was appropriate.
	 */
	public boolean setGraphScale(double scale) {
		double newScale = scale;
		if (scale > ZOOM_IN_LIMIT) {
			newScale = ZOOM_IN_LIMIT;
			zoomInButton.setEnabled(false);
		} else {
			zoomInButton.setEnabled(true);
		}

		if (scale < ZOOM_OUT_LIMIT) {
			newScale = ZOOM_OUT_LIMIT;
			zoomOutButton.setEnabled(false);
		} else {
			zoomOutButton.setEnabled(true);
		}

		graph.setScale(newScale);
		graphScroller.repaint();
		return newScale == scale;
	}

	private void showPopupMenu(MouseEvent e) {
		// right click
		List<AbstractGraphVertex> chosenCells = graph.getVertexesAtPoint(e.getPoint());
		if (chosenCells != null && chosenCells.size() > 0) {

			List<DataItem> items = new ArrayList<DataItem>();
			for (Dataset bean : application.getSelectionManager().getSelectedDataBeans()) {
				items.add(bean);
			}
			application.showPopupMenuFor(e, items);

		} else {
			Dataset nullBean = null;
			application.showPopupMenuFor(e, nullBean);
		}
	}

	// used to send values for the new thread in method zoomInToolUsed
	private Point center;

	/**
	 * Zoom in
	 * 
	 * @param e
	 */
	private void zoomInToolUsed() {
		// If the zoom-limit has been reached, don't do anything
		if (this.setGraphScale(graph.getScale() * this.ZOOM_FACTOR)) {

			JViewport view = this.getScroller().getViewport();
			Point viewPos = view.getViewPosition();
			center = new Point((int) (viewPos.getX() + view.getWidth() / 2), (int) (viewPos.getY() + view.getHeight() / 2));

			logger.debug("view.getX: " + view.getX() + ", viewWidth: " + view.getWidth());

			// Make the coordinates to correspond the new size of the graph
			center.x = (int) (center.getX() * this.ZOOM_FACTOR);
			center.y = (int) (center.getY() * this.ZOOM_FACTOR);

			/*
			 * The new scale won't be applied before this method is ended. The new scale changes the size of the graph, which affects on the
			 * calculations done in method pointToCenter. We couuld estimate the new size of the graph, but that doesn't help because we
			 * still cant scroll to the new areas.
			 */

			Runnable centerer = new Thread() {
				public void run() {
					try {
						Thread.sleep(10);
					} catch (InterruptedException e) {
						// Nothing serious happened
					}
					GraphPanel.this.pointToCenter(GraphPanel.this.center);
				}
			};

			updateGridSize();

			// try it to avoid flickering when possible
			this.pointToCenter(GraphPanel.this.center);

			// and make sure it's done in every situation;
			SwingUtilities.invokeLater(centerer);
		}
		graph.repaint();
	}

	/**
	 * Zooms out.
	 */
	private void zoomOutToolUsed() {
		if (this.setGraphScale(graph.getScale() / this.ZOOM_FACTOR)) {
			JViewport view = this.getScroller().getViewport();
			Point viewPos = view.getViewPosition();
			center = new Point((int) (viewPos.getX() + view.getWidth() / 2), (int) (viewPos.getY() + view.getHeight() / 2));

			center.x = (int) (center.getX() / this.ZOOM_FACTOR);
			center.y = (int) (center.getY() / this.ZOOM_FACTOR);

			updateGridSize();

			this.pointToCenter(center);
		}

		this.repaint();
		this.getScroller().repaint();
		graph.repaint();
	}

	private void updateGridSize() {
		double size = DEFAULT_GRID_SIZE;
		if (graph.getScale() <= 0.5) {
			size = BIGGER_GRID_SIZE;
		}
		if (graph.getScale() <= 0.25) {
			size = HUGE_GRID_SIZE;
		}

		graph.setGridSize(size);
	}

	/**
	 * Enables the history button if only one dataset is selected
	 * 
	 * @see java.beans.PropertyChangeListener#propertyChange(java.beans.PropertyChangeEvent)
	 */
	public void propertyChange(PropertyChangeEvent e) {

		if (e instanceof DatasetChoiceEvent) {

			SwingTools.runInEventDispatchThread(new Runnable() {
				public void run() {
					if (application.getSelectionManager().getSelectedDataBean() != null) {
						Dataset bean = application.getSelectionManager().getSelectedDataBean();
						graph.scrollCellToVisibleAnimated(graph.getVertexMap().get(bean));
					}

					graph.repaint();
					historyButton.setEnabled(application.getSelectionManager().getSelectedItem() instanceof Dataset);
				}
			});

		}
	}

	public void dataChanged(DataChangeEvent event) {
		if (!(event instanceof ContentChangedEvent || event instanceof LinksChangedEvent || !(event instanceof DataItemCreatedEvent))) {
			SwingTools.runInEventDispatchThread(new Runnable() {
				public void run() {
					autoZoom();
				}
			});
		}
	}

	public class WorkflowSelectionListener implements GraphSelectionListener {

		public void valueChanged(GraphSelectionEvent e) {
			if (!internalSelection) {

				boolean emptySelection = (graph.getSelectionCount() == 0);

				application.getSelectionManager().clearAll(emptySelection, graph);

				ArrayList<AbstractGraphVertex> vertexes = new ArrayList<AbstractGraphVertex>();

				for (Object obj : graph.getSelectionCells()) {

					if (obj instanceof GroupVertex) {
						GroupVertex group = (GroupVertex) obj;

						vertexes.addAll(group.getChildVertexes());
					}

					if (obj instanceof GraphVertex) {
						vertexes.add((AbstractGraphVertex) obj);
					}
				}

				ArrayList<DataItem> items = new ArrayList<DataItem>();

				for (AbstractGraphVertex vertex : vertexes) {
					items.add((DataItem) vertex.getUserObject());
				}

				application.getSelectionManager().selectMultiple(items, graph);
			}
		}
	}
}
