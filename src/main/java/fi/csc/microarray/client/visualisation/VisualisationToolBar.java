package fi.csc.microarray.client.visualisation;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.ToolBarComponentFactory;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * This panel contains the options for different data visualizations. It is a
 * lower left part of the operations panel, and controls what is shown on the
 * results panel.
 * 
 * @author Janne Käki, akallio
 * 
 */
public class VisualisationToolBar extends JToolBar implements ActionListener, PropertyChangeListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(VisualisationToolBar.class);

	public boolean isMaximised = false;

	private ClientApplication application = Session.getSession().getApplication();

	private JButton redrawButton = ToolBarComponentFactory.createButton("Redraw", VisualConstants.REDRAW_ICON, true, false);
	private JButton helpButton = ToolBarComponentFactory.createButton("Help", VisualConstants.QUESTION_MARK_ICON, true, false);
	private JButton maximiseButton = ToolBarComponentFactory.createButton("Maximise", VisualConstants.MAXIMISE_ICON, true, false);
	private JButton splitButton = ToolBarComponentFactory.createButton("Duplicate", VisualConstants.SPLIT_ICON, true, false);
	private JButton detachButton = ToolBarComponentFactory.createButton("Detach", VisualConstants.TO_WINDOW_ICON, true, false);

	// private VisualisationListModel methodListModel = new
	// VisualisationListModel();
	private JComboBox methodChoiceBox = ToolBarComponentFactory.createComboBox();

	private String helpAddress;

	private boolean userComboAction;

	private boolean isSplit;

	public VisualisationToolBar() {

		super();

		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		this.setFloatable(false);
		this.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);

		Dimension ZERO_SIZE = new Dimension(0, 0);

		redrawButton.setMinimumSize((Dimension) ZERO_SIZE.clone());
		helpButton.setMinimumSize((Dimension) ZERO_SIZE.clone());
		maximiseButton.setMinimumSize((Dimension) ZERO_SIZE.clone());
		splitButton.setMinimumSize((Dimension) ZERO_SIZE.clone());
		detachButton.setMinimumSize((Dimension) ZERO_SIZE.clone());

		methodChoiceBox.setBackground(Color.WHITE);
		methodChoiceBox.addActionListener(this);
		// methodChoiceBox.setPreferredSize(ToolBarComponentFactory.COMBOBOX_SIZE);
		methodChoiceBox.setRenderer(new ComboBoxRenderer());

		JPanel buttonPanel = new JPanel(new GridLayout(1, 4));
		buttonPanel.setOpaque(false);
		buttonPanel.add(helpButton);
		buttonPanel.add(redrawButton);
		buttonPanel.add(maximiseButton);
		// buttonPanel.add(splitButton); // splitting disabled, needs more
		// thinking
		buttonPanel.add(detachButton);

		JLabel titleLabel = new JLabel("Method: ");

		this.add(titleLabel);
		this.add(methodChoiceBox);
		this.add(Box.createHorizontalStrut(3));
		this.add(Box.createHorizontalGlue());
		this.add(buttonPanel);

		redrawButton.setEnabled(false);
		redrawButton.addActionListener(this);
		helpButton.addActionListener(this);
		maximiseButton.addActionListener(this);
		splitButton.addActionListener(this);
		detachButton.addActionListener(this);

		helpButton.setVisible(false);

		refreshVisualisationList(VisualisationMethod.NONE, null);

		// start listening
		application.addPropertyChangeListener(this);
	}

	public Vector<Component> getFocusComponents() {
		Vector<Component> order = new Vector<Component>();
		order.add(methodChoiceBox);
		return order;
	}

	/**
	 * Can be called also outside EDT
	 * 
	 * @param helpAddress
	 */
	public void setHelpAddress(String helpAddress) {

		this.helpAddress = helpAddress;
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				helpButton.setVisible(VisualisationToolBar.this.helpAddress != null);
				VisualisationToolBar.this.repaint();
			}
		});
	}

	private void refreshVisualisationList(VisualisationMethod method, List<DataBean> datas) {

		// update maximise button
		maximiseButton.setEnabled(datas != null && datas.size() > 0);

		// update method list
		fillMethodsFor(datas);
		methodChoiceBox.setEnabled(datas != null && datas.size() > 0);
		redrawButton.setEnabled(method != VisualisationMethod.NONE);
		splitButton.setEnabled(method != VisualisationMethod.NONE || isSplit);
		detachButton.setEnabled(method != VisualisationMethod.NONE);
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to
	 * listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		VisualisationFrame visualisation;

		// TODO Find a reason for the stack overflow and remove this if and else
		// block
		if (application.getVisualisationFrameManager() != null) {
			visualisation = application.getVisualisationFrameManager().getFrame(FrameType.MAIN);
		} else {
			// Only initialising, no real actions
			return;
		}

		if (source == methodChoiceBox) {
			if (userComboAction) {

				VisualisationMethod method = (VisualisationMethod) methodChoiceBox.getSelectedItem();
				if (method == null) {
					method = VisualisationMethod.NONE;
				}

				if (visualisation.getMethod() != method) {
					application.setVisualisationMethod(method, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
				}
				redrawButton.setEnabled(method != VisualisationMethod.NONE);
			}

		} else if (source == maximiseButton) {
			maximiseOrRestoreVisualisation();

		} else if (source == redrawButton) {
			application.setVisualisationMethod(visualisation.getMethod(), visualisation.getVariables(), visualisation.getDatas(), visualisation.getType());
			// FIXME redraw second visualisation if present
		} else if (source == helpButton) {
			viewHelp();
		} else if (source == splitButton) {
			split();

		} else if (source == detachButton) {
			detach();
		}
	}

	public void maximiseOrRestoreVisualisation() {
		application.setMaximisedVisualisationMode(isMaximised = !isMaximised); // flips
																				// isMaximised
																				// and
																				// uses
																				// the
																				// new
																				// value
		maximiseButton.setText(getMaximiseButtonText());
		maximiseButton.setIcon(getMaximiseButtonIcon());
	}

	public Icon getMaximiseButtonIcon() {
		return isMaximised ? VisualConstants.RESTORE_ICON : VisualConstants.MAXIMISE_ICON;
	}

	public String getMaximiseButtonText() {
		return isMaximised ? "Restore" : "Maximise";
	}

	public void viewHelp() {
		application.viewHelp(helpAddress);
	}

	public void detach() {
		VisualisationFrame visualisation = application.getVisualisationFrameManager().getFrame(FrameType.MAIN);

		application.setVisualisationMethod(visualisation.getMethod(), visualisation.getVariables(), visualisation.getDatas(), FrameType.WINDOW);
	}

	public void split() {

		VisualisationFrame visualisation = application.getVisualisationFrameManager().getFrame(FrameType.MAIN);

		if (isSplit = !isSplit) {
			application.setVisualisationMethod(visualisation.getMethod(), visualisation.getVariables(), visualisation.getDatas(), FrameType.SIDE);

		} else {
			application.getVisualisationFrameManager().closeAllByType(FrameType.SIDE);

			splitButton.setEnabled(!methodChoiceBox.getSelectedItem().equals(VisualisationMethod.NONE));
		}

		splitButton.setText(getSplitText());
		splitButton.setIcon(getSplitIcon());
	}

	public String getSplitText() {
		return isSplit ? "Close second" : "Duplicate";
	}

	public Icon getSplitIcon() {
		return isSplit ? VisualConstants.CLOSE_ICON : VisualConstants.SPLIT_ICON;
	}

	public void propertyChange(PropertyChangeEvent event) {

		if (event instanceof VisualisationMethodChangedEvent) {


			VisualisationMethodChangedEvent e = (VisualisationMethodChangedEvent) event;

			// FIXME multiple window compatibility
			if (e.getTarget() == FrameType.MAIN) {
				
				// update help button
				setHelpAddress(e.getNewMethod().getHelpAddress());
				
				refreshVisualisationList(e.getNewMethod(), e.getDatas());
				try {

					if (e.getNewMethod() != null) {
						logger.debug("updating GUI to method " + e.getNewMethod());

						userComboAction = false;
						methodChoiceBox.setSelectedItem(e.getNewMethod());
						userComboAction = true;

						// Make sure the selected visualisation is shown
						methodChoiceBox.repaint();
					}
				} catch (Exception exc) {
					application.reportException(exc);
				}
			}
		} else if (event instanceof DatasetChoiceEvent) {
			List<DataBean> currentDatas = application.getSelectionManager().getSelectedDataBeans();
			List<DataBean> newDatas = application.getVisualisationFrameManager().getFrame(FrameType.MAIN).getDatas();

			// If same
			if (currentDatas == null || newDatas == null || !(currentDatas.containsAll(newDatas) && newDatas.containsAll(currentDatas))) {

				application.setVisualisationMethod(VisualisationMethod.NONE, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);

				refreshVisualisationList(VisualisationMethod.NONE, application.getSelectionManager().getSelectedDataBeans());
			}
		}
	}

	public void fillMethodsFor(List<DataBean> datas) {
		// Arrays.asList doesn't support removing, so we need a new one
		List<VisualisationMethod> applicableVisualisations = new ArrayList<VisualisationMethod>();
		applicableVisualisations.addAll(Session.getSession().getVisualisations().getVisualisationMethods());

		List<VisualisationMethod> onlyNoneList = new ArrayList<VisualisationMethod>();
		onlyNoneList.add(VisualisationMethod.NONE);

		if (datas != null) {

			for (VisualisationMethod method : Session.getSession().getVisualisations().getVisualisationMethods()) {

				try {

					if (datas.size() == 1) {
						if (!method.isApplicableTo(datas.get(0))) {
							applicableVisualisations.remove(method);
						}
					} else {

						if (!method.isApplicableTo(datas)) {
							applicableVisualisations.remove(method);
						}
					}

				} catch (Exception e) {
					application.reportException(new MicroarrayException("Unable to check applicable visualisations for the dataset", e));

					applicableVisualisations = onlyNoneList;
					break;
				}
			}

			userComboAction = false;
			methodChoiceBox.removeAllItems();
			Visualisation.fillCompoBox(methodChoiceBox, applicableVisualisations.toArray());
			userComboAction = true;

		} else {
			applicableVisualisations = onlyNoneList;

			userComboAction = false;
			methodChoiceBox.removeAllItems();
			Visualisation.fillCompoBox(methodChoiceBox, applicableVisualisations.toArray());
			userComboAction = true;
		}
	}
}
