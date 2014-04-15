package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;

import net.miginfocom.swing.MigLayout;

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
import fi.csc.microarray.module.basic.BasicModule.VisualisationMethods;

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

	private JButton helpButton = ToolBarComponentFactory.createButton("Help", VisualConstants.QUESTION_MARK_ICON, true, false);
	private JButton maximiseButton = ToolBarComponentFactory.createButton("Maximise", VisualConstants.MAXIMISE_ICON, true, false);
	private JButton detachButton = ToolBarComponentFactory.createButton("Detach", VisualConstants.TO_WINDOW_ICON, true, false);
	private JButton closeButton = ToolBarComponentFactory.createButton("Close", VisualConstants.CLOSE_ICON, true, false);

	JPanel buttonPanel;

	private JComboBox<VisualisationMethod> methodChoiceBox = ToolBarComponentFactory.createComboBox();

	private String helpAddress;

	private boolean userComboAction;


	public VisualisationToolBar() {

		super();

		this.setLayout(new MigLayout("alignx right, height 22!, insets 0"));
		this.setFloatable(false);
		this.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);

		methodChoiceBox.addActionListener(this);
		methodChoiceBox.setRenderer(new ComboBoxRenderer());

		String width = "width 100";
		
		this.add(methodChoiceBox, "width 200, pushx, aligny top");
		this.add(helpButton, width);
		this.add(maximiseButton, width);
		this.add(detachButton, width);
		this.add(closeButton, width);


		helpButton.addActionListener(this);
		maximiseButton.addActionListener(this);
		detachButton.addActionListener(this);
		closeButton.addActionListener(this);

		helpButton.setVisible(false);

		refreshVisualisationList(VisualisationMethod.getDefault(), null);

		// start listening
		application.addClientEventListener(this);
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
		if (method != VisualisationMethods.DATA_DETAILS) {
			methodChoiceBox.setVisible(true);
			fillMethodsFor(datas);
			methodChoiceBox.setEnabled(datas != null && datas.size() > 0);
		} else {
			methodChoiceBox.setVisible(false);
		}
			
		closeButton.setEnabled(method != VisualisationMethods.DATA_DETAILS);
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
					method = VisualisationMethod.getDefault();
				}

				if (visualisation.getMethod() != method) {
					application.setVisualisationMethod(method, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
				}
			}

		} else 
			if (source == maximiseButton) {
			maximiseOrRestoreVisualisation();
		} else if (source == helpButton) {
			viewHelp();
		} else if (source == detachButton) {
			detach();
		} else if (source == closeButton) {
			application.setVisualisationMethod(VisualisationMethods.DATA_DETAILS, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
		}
	}

	public void maximiseOrRestoreVisualisation() {
		application.setMaximisedVisualisationMode(isMaximised = !isMaximised); // flips
		// isMaximised and uses the new value
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

			VisualisationMethod method = null;
			method = VisualisationMethod.getDefault();			
			
			if (currentDatas == null || newDatas == null || !(currentDatas.containsAll(newDatas) && newDatas.containsAll(currentDatas))) {

				application.setVisualisationMethod(method, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);

				refreshVisualisationList(method, application.getSelectionManager().getSelectedDataBeans());
			}
		}
	}

	public static List<VisualisationMethod> getMethodsFor(List<DataBean> datas) {
		// Arrays.asList doesn't support removing, so we need a new one
		List<VisualisationMethod> applicableVisualisations = new ArrayList<VisualisationMethod>();
		applicableVisualisations.addAll(Session.getSession().getVisualisations().getVisualisationMethods());

		List<VisualisationMethod> onlyDefaultList = new ArrayList<VisualisationMethod>();
		onlyDefaultList.add(VisualisationMethod.getDefault());

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
					Session.getSession().getApplication().reportException(new MicroarrayException("Unable to check applicable visualisations for the dataset", e));

					applicableVisualisations = onlyDefaultList;
					break;
				}
			}
		} else {
			applicableVisualisations = onlyDefaultList;
		}

		return applicableVisualisations;
	}

	public void fillMethodsFor(List<DataBean> datas) {

		userComboAction = false;
		methodChoiceBox.removeAllItems();
		Visualisation.fillComboBox(methodChoiceBox, getMethodsFor(datas).toArray(new VisualisationMethod[0]));
		userComboAction = true;
	}
}