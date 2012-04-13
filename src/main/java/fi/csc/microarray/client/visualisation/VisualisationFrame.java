package fi.csc.microarray.client.visualisation;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.SwingConstants;

import org.testng.log4testng.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataItemRemovedEvent;

public abstract class VisualisationFrame implements DataChangeListener {

	public abstract void setContent(JComponent visualisationComponent);

	private static final Color BG = Color.white;

	protected ClientApplication application = Session.getSession().getApplication();

	private static final String WAIT_PANEL_NAME = "wait";
	private static final String VISUALISATION_PANEL_NAME = "visualisation";

	private CardLayout viewChangerLayout = new CardLayout();
	private JPanel viewChangerPanel = new JPanel(viewChangerLayout);
	private JComponent lastVisualisationPanel;
	private JSplitPane paramSplit;

	private VisualisationMethod method;

	private List<Variable> variables;

	private List<DataBean> datas;

	protected FrameType type;

	private Visualisation visualiser;

	private static final Logger logger = Logger.getLogger(VisualisationFrame.class);

	/**
	 * Creates a new VisualisationFrame with a CardLayout-powered view changer.
	 * 
	 * @throws IOException
	 */
	public VisualisationFrame(FrameType type) {
		this.type = type;

		// initialise panels
		viewChangerPanel.setBackground(BG);

		// initialise wait panel
		JPanel waitPanel = new JPanel(new BorderLayout());
		JLabel waitLabel = new JLabel("Visualising, please wait...");
		waitLabel.setFont(waitLabel.getFont().deriveFont(Font.BOLD));
		waitLabel.setHorizontalAlignment(SwingConstants.CENTER);
		waitPanel.add(waitLabel, BorderLayout.CENTER);
		viewChangerPanel.add(waitPanel, WAIT_PANEL_NAME);

		// show empty panel
		try {
			this.showVisualisationComponent(VisualisationMethod.NONE.getVisualiser(this).getVisualisation(new ArrayList<DataBean>()));
		} catch (Exception e) {
			application.reportException(e);
		}
		this.setContent(viewChangerPanel);
		viewChangerLayout.show(viewChangerPanel, VISUALISATION_PANEL_NAME);

		application.getDataManager().addDataChangeListener(this);
	}

	private class SplitSizeHandler implements PropertyChangeListener {
		public void propertyChange(PropertyChangeEvent e) {
			int leftLimit = (int) paramSplit.getWidth() - (int) Visualisation.PARAMETER_SIZE.getWidth();
			if (paramSplit.getDividerLocation() < leftLimit) {
				paramSplit.setDividerLocation(leftLimit);
			}
		}
	}

	public JComponent createVisualisation(VisualisationMethodChangedEvent e) {

		// Create new visualiser only if needed to keep the settings made in settings panel
		if (this.datas != e.getDatas() || this.method != e.getNewMethod()) {
			this.datas = e.getDatas();
			this.method = e.getNewMethod();

			removeVisualiser();
			visualiser = method.getVisualiser(this);
		}
		this.variables = e.getVariables();

		// parameter panel has to be first one to make it initialised before the
		// data is set (scatterplot)
		JPanel parametersPanel = visualiser.getParameterPanel();
		logger.debug("parametersPanel for method " + method + " contains: " + parametersPanel);
		if (parametersPanel != null) {
			paramSplit = new JSplitPane();
			parametersPanel.setMinimumSize(new Dimension(0, 0));
			paramSplit.setRightComponent(parametersPanel);
			// To show the width limit of parameter panel
			paramSplit.setContinuousLayout(true);
			// To keep the parameter panel size constant
			paramSplit.setResizeWeight(1.0);

			SplitSizeHandler sizeHandler = new SplitSizeHandler();
			paramSplit.addPropertyChangeListener(JSplitPane.DIVIDER_LOCATION_PROPERTY, sizeHandler);
		}

		JComponent visualisationComponent = null;

		JComponent componentToReturn = null;
		try {
			if (visualiser.isForMultipleDatas()) {
				visualisationComponent = visualiser.getVisualisation(datas);
			} else if (visualiser.isForSingleData()) {
				DataBean data = datas.size() > 0 ? datas.get(0) : null;
				visualisationComponent = visualiser.getVisualisation(data);
			}

			if (parametersPanel != null) {
				paramSplit.setLeftComponent(visualisationComponent);
				componentToReturn = paramSplit;
			} else {
				componentToReturn = visualisationComponent;
			}

			if (logger.isDebugEnabled()) {
				logger.debug("visualisationpanel contains following components:");
				for (Component component : paramSplit.getComponents()) {
					logger.debug("\t" + component);
				}
			}
		} catch (Exception e1) {
			application.reportException(e1);
			componentToReturn = visualiser.getDefaultVisualisation();
		}

		return componentToReturn;
	}

	/**
	 * Is not explicitly synchronised, should be called only from EDT.
	 */
	public void showVisualisationComponent(JComponent panel) {
		// Clean
		removeVisualisationComponent();

		if (panel == null) {
			try {
				panel = VisualisationMethod.NONE.getVisualiser(this).getVisualisation(new ArrayList<DataBean>());
			} catch (Exception e) {
				application.reportException(e);
			}
		}

		this.lastVisualisationPanel = panel;

		String title = "Visualisation";
		if (datas != null && datas.size() > 0) {
			title += " of ";
			for (int i = 0; i < datas.size(); i++) {
				title += datas.get(i);
				if (i != datas.size() - 1) {
					title += (i == datas.size() - 2) ? " and " : ", ";
				}
			}
		}

		this.setTitle(title);

		// Do
		viewChangerPanel.add(panel, VISUALISATION_PANEL_NAME);
		viewChangerLayout.show(viewChangerPanel, VISUALISATION_PANEL_NAME);
		// Split obeys divider locations only after it's shown, else side visualisations hide parameters
		if (paramSplit != null) {
			paramSplit.setDividerLocation(0.5);
		}
	}

	/**
	 * Is not explicitly synchronised, should be called only from EDT.
	 */
	public void showWaitPanel() {
		removeVisualisationComponent();
		viewChangerLayout.show(viewChangerPanel, WAIT_PANEL_NAME);
	}

	/**
	 * Removing visualisation panels immediately is crucial for memory efficiency. Is not explicitly synchronised, should be called only
	 * from EDT.
	 */
	public void removeVisualisationComponent() {

		if (lastVisualisationPanel != null) {
			
			// remove all references to visualisation panel
			viewChangerPanel.remove(lastVisualisationPanel);
			viewChangerLayout.removeLayoutComponent(lastVisualisationPanel);

			lastVisualisationPanel = null;
		}
	}

	private void removeVisualiser() {
		if (visualiser instanceof PropertyChangeListener) {
			application.removeClientEventListener((PropertyChangeListener) visualiser);
		}
		if (visualiser != null) {
			visualiser.removeVisualisation();
		}
	}

	public VisualisationMethod getVisualisationMethod() {
		return method;
	}

	public List<Variable> getVariables() {
		return variables;
	}

	public List<DataBean> getDatas() {
		return datas;
	}

	public FrameType getType() {
		return type;
	}

	public VisualisationMethod getMethod() {
		return method;
	}

	public void dataChanged(DataChangeEvent e) {
		
		if (e instanceof DataItemRemovedEvent) {
			// Data removed, check if it was part of visualisation
			if (datas != null && datas.contains(((DataItemRemovedEvent) e).getDataItem())) {
				application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, VisualisationMethod.NONE, null, null, type, this));
			}
			
			// Data removed, check if context links should be refreshed
			if (method == null || method == VisualisationMethod.NONE) {
				updateContextLinks();
			}

			
		} else if (e instanceof DataItemCreatedEvent) {

			// Data added, check if context links should be refreshed
			if (method == null || method == VisualisationMethod.NONE) {
				updateContextLinks();
			}
		}
	}

	protected abstract void updateContextLinks();
	
	void setTitle(String title) {
		return;
	}
}
