package com.nexes.wizard;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Insets;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.HashMap;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

import fi.csc.microarray.client.VisualConstants;

/**
 * This class implements a basic wizard dialog, where the programmer can
 * insert one or more Components to act as panels. These panels can be navigated
 * through arbitrarily using the 'Next' or 'Back' buttons, or the dialog itself
 * can be closed using the 'Cancel' button. Note that even though the dialog
 * uses a CardLayout manager, the order of the panels is not linear. Each panel
 * determines at runtime what its next and previous panel will be.
 */
public class Wizard extends WindowAdapter implements PropertyChangeListener {

	/**
	 * Indicates that the 'Finish' button was pressed to close the dialog.
	 */
	public static final int FINISH_RETURN_CODE = 0;

	/**
	 * Indicates that the 'Cancel' button was pressed to close the dialog, or
	 * the user pressed the close box in the corner of the window.
	 */
	public static final int CANCEL_RETURN_CODE = 1;

	/**
	 * Indicates that the dialog closed due to an internal error.
	 */
	public static final int ERROR_RETURN_CODE = 2;

	/**
	 * The String-based action command for the 'Next' button.
	 */
	public static final String NEXT_BUTTON_ACTION_COMMAND = "NextButtonActionCommand";

	/**
	 * The String-based action command for the 'Back' button.
	 */
	public static final String BACK_BUTTON_ACTION_COMMAND = "BackButtonActionCommand";

	/**
	 * The String-based action command for the 'Cancel' button.
	 */
	public static final String CANCEL_BUTTON_ACTION_COMMAND = "CancelButtonActionCommand";

	/**
	 * The String-based action command for the 'Finish' button.
	 */
	public static final String FINISH_BUTTON_ACTION_COMMAND = "FinishButtonActionCommand";

	public static final String START_NAVIGATION_PANEL = "StartNavigationPanel";

	public static final String FINISH_NAVIGATION_PANEL = "FinishNavigationPanel";

	public static final String NEXT_NAVIGATION_PANEL = "NextNavigationPanel";

	public static final String FINISH_NEXT_NAVIGATION_PANEL = "FinishNextNavigationPanel";

	// The i18n text used for the buttons. Loaded from a property resource file.    

	static String BACK_TEXT = "Back";

	static String NEXT_TEXT = "Next";

	static String FINISH_TEXT = "Finish";

	static String CANCEL_TEXT = "Cancel";

	private WizardModel wizardModel;

	private WizardController wizardController;

	private JDialog wizardDialog;

	private JPanel contentPanel;

	private CardLayout contentLayout;

	private JPanel buttonPanel;

	private CardLayout buttonLayout;

	private Map<String, JButton> backButtons = new HashMap<String, JButton>();

	private Map<String, JButton> nextButtons = new HashMap<String, JButton>();

	private Map<String, JButton> cancelButtons = new HashMap<String, JButton>();

	private Map<String, JButton> finishButtons = new HashMap<String, JButton>();

	private int returnCode;

	/**
	 * Default constructor. This method creates a new WizardModel object and passes it
	 * into the overloaded constructor.
	 */
	public Wizard() {
		this((Frame) null);
	}

	/**
	 * This method accepts a java.awt.Dialog object as the javax.swing.JDialog's
	 * parent.
	 * @param owner The java.awt.Dialog object that is the owner of this dialog.
	 */
	public Wizard(Dialog owner) {
		wizardModel = new WizardModel();
		wizardDialog = new JDialog(owner);
		initComponents();
	}

	/**
	 * This method accepts a java.awt.Frame object as the javax.swing.JDialog's
	 * parent.
	 * @param owner The java.awt.Frame object that is the owner of the javax.swing.JDialog.
	 */
	public Wizard(Frame owner) {
		wizardModel = new WizardModel();
		wizardDialog = new JDialog(owner);
		initComponents();
		if (owner != null) {
			wizardDialog.setLocationRelativeTo(owner);
		}
	}

	/**
	 * Returns an instance of the JDialog that this class created. This is useful in
	 * the event that you want to change any of the JDialog parameters manually.
	 * @return The JDialog instance that this class created.
	 */
	public JDialog getDialog() {
		return wizardDialog;
	}

	/**
	 * Returns the owner of the generated javax.swing.JDialog.
	 * @return The owner (java.awt.Frame or java.awt.Dialog) of the javax.swing.JDialog generated
	 * by this class.
	 */
	public Component getOwner() {
		return wizardDialog.getOwner();
	}

	/**
	 * Sets the title of the generated javax.swing.JDialog.
	 * @param s The title of the dialog.
	 */
	public void setTitle(String s) {
		wizardDialog.setTitle(s);
	}

	/**
	 * Returns the current title of the generated dialog.
	 * @return The String-based title of the generated dialog.
	 */
	public String getTitle() {
		return wizardDialog.getTitle();
	}

	/**
	 * Sets the modality of the generated javax.swing.JDialog.
	 * @param b the modality of the dialog
	 */
	public void setModal(boolean b) {
		wizardDialog.setModal(b);
	}

	/**
	 * Returns the modality of the dialog.
	 * @return A boolean indicating whether or not the generated javax.swing.JDialog is modal.
	 */
	public boolean isModal() {
		return wizardDialog.isModal();
	}

	/**
	 * Convienence method that displays a modal wizard dialog and blocks until the dialog
	 * has completed.
	 * @return Indicates how the dialog was closed. Compare this value against the RETURN_CODE
	 * constants at the beginning of the class.
	 */
	public int showModalDialog() {

		wizardDialog.setModal(true);
		wizardDialog.pack();
		wizardDialog.setVisible(true);

		return returnCode;
	}

	/**
	 * Returns the current model of the wizard dialog.
	 * @return A WizardModel instance, which serves as the model for the wizard dialog.
	 */
	public WizardModel getModel() {
		return wizardModel;
	}

	/**
	 * Add a Component as a panel for the wizard dialog by registering its
	 * WizardPanelDescriptor object. Each panel is identified by a unique Object-based
	 * identifier (often a String), which can be used by the setCurrentPanel()
	 * method to display the panel at runtime.
	 * @param id An Object-based identifier used to identify the WizardPanelDescriptor object.
	 * @param panel The WizardPanelDescriptor object which contains helpful information about the panel.
	 */
	public void registerWizardPanel(Object id, WizardPanelDescriptor panel) {

		//  Add the incoming panel to our JPanel display that is managed by
		//  the CardLayout layout manager.
		//  Add the incoming panel to our JPanel display that is managed by
		//  the CardLayout layout manager.

		contentPanel.add(panel.getPanelComponent(), id);

		//  Set a callback to the current wizard.

		panel.setWizard(this);

		//  Place a reference to it in the model. 

		wizardModel.registerPanel(id, panel);
	}

	/**
	 * Displays the panel identified by the object passed in. This is the same Object-based
	 * identified used when registering the panel.
	 * @param id The Object-based identifier of the panel to be displayed.
	 */
	public void setCurrentPanel(Object id) {

		//  Get the hashtable reference to the panel that should
		//  be displayed. If the identifier passed in is null, then close
		//  the dialog.

		if (id == null)
			close(ERROR_RETURN_CODE);

		WizardPanelDescriptor oldPanelDescriptor = wizardModel
				.getCurrentPanelDescriptor();
		if (oldPanelDescriptor != null)
			oldPanelDescriptor.aboutToHidePanel();

		wizardModel.setCurrentPanel(id);
		wizardModel.getCurrentPanelDescriptor().aboutToDisplayPanel();

		//  Show the content panel in the dialog.

		contentLayout.show(contentPanel, id.toString());

		// Select and show the correct button panel

		String panelType = wizardModel.getCurrentPanelDescriptor()
				.getButtonPanelType();
		buttonLayout.show(buttonPanel, panelType);

		wizardModel.getCurrentPanelDescriptor().displayingPanel();
	}

	/**
	 * Method used to listen for property change events from the model and update the
	 * dialog's graphical components as necessary.
	 * @param evt PropertyChangeEvent passed from the model to signal that one of its properties has changed value.
	 */
	public void propertyChange(PropertyChangeEvent evt) {
		String panelType = wizardModel.getCurrentPanelDescriptor()
				.getButtonPanelType();
		JButton button;

		if (evt.getPropertyName().equals(
				WizardModel.CURRENT_PANEL_DESCRIPTOR_PROPERTY)) {
			wizardController.resetButtonsToPanelRules();
		} else if (evt.getPropertyName().equals(
				WizardModel.NEXT_BUTTON_TEXT_PROPERTY)) {
			button = (JButton) nextButtons.get(panelType);
			if (button != null)
				button.setText(evt.getNewValue().toString());
		} else if (evt.getPropertyName().equals(
				WizardModel.FINISH_BUTTON_TEXT_PROPERTY)) {
			button = (JButton) finishButtons.get(panelType);
			if (button != null)
				button.setText(evt.getNewValue().toString());
		} else if (evt.getPropertyName().equals(
				WizardModel.BACK_BUTTON_TEXT_PROPERTY)) {
			button = (JButton) backButtons.get(panelType);
			if (button != null)
				button.setText(evt.getNewValue().toString());
		} else if (evt.getPropertyName().equals(
				WizardModel.CANCEL_BUTTON_TEXT_PROPERTY)) {
			button = (JButton) cancelButtons.get(panelType);
			if (button != null)
				button.setText(evt.getNewValue().toString());
		} else if (evt.getPropertyName().equals(
				WizardModel.NEXT_BUTTON_ENABLED_PROPERTY)) {
			button = (JButton) nextButtons.get(panelType);
			if (button != null)
				button.setEnabled(((Boolean) evt.getNewValue()).booleanValue());
		} else if (evt.getPropertyName().equals(
				WizardModel.FINISH_BUTTON_ENABLED_PROPERTY)) {
			button = (JButton) finishButtons.get(panelType);
			if (button != null)
				button.setEnabled(((Boolean) evt.getNewValue()).booleanValue());
		} else if (evt.getPropertyName().equals(
				WizardModel.BACK_BUTTON_ENABLED_PROPERTY)) {
			button = (JButton) backButtons.get(panelType);
			if (button != null)
				button.setEnabled(((Boolean) evt.getNewValue()).booleanValue());
		} else if (evt.getPropertyName().equals(
				WizardModel.CANCEL_BUTTON_ENABLED_PROPERTY)) {
			button = (JButton) cancelButtons.get(panelType);
			if (button != null)
				button.setEnabled(((Boolean) evt.getNewValue()).booleanValue());
		}
	}

	/**
	 * Retrieves the last return code set by the dialog.
	 * @return An integer that identifies how the dialog was closed. See the *_RETURN_CODE
	 * constants of this class for possible values.
	 */
	public int getReturnCode() {
		return returnCode;
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @return A boolean indicating if the button is enabled.
	 */
	public boolean getBackButtonEnabled() {
		return wizardModel.getBackButtonEnabled().booleanValue();
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @param boolean newValue The new enabled status of the button.
	 */
	public void setBackButtonEnabled(boolean newValue) {
		wizardModel.setBackButtonEnabled(new Boolean(newValue));
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @return A boolean indicating if the button is enabled.
	 */
	public boolean getNextButtonEnabled() {
		return wizardModel.getNextButtonEnabled().booleanValue();
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @param boolean newValue The new enabled status of the button.
	 */
	public void setNextButtonEnabled(boolean newValue) {
		wizardModel.setNextButtonEnabled(new Boolean(newValue));
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @return A boolean indicating if the button is enabled.
	 */
	public boolean getFinishButtonEnabled() {
		return wizardModel.getFinishButtonEnabled().booleanValue();
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @param boolean newValue The new enabled status of the button.
	 */
	public void setFinishButtonEnabled(boolean newValue) {
		wizardModel.setFinishButtonEnabled(new Boolean(newValue));
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @return A boolean indicating if the button is enabled.
	 */
	public boolean getCancelButtonEnabled() {
		return wizardModel.getCancelButtonEnabled().booleanValue();
	}

	/**
	 * Mirrors the WizardModel method of the same name.
	 * @param boolean newValue The new enabled status of the button.
	 */
	public void setCancelButtonEnabled(boolean newValue) {
		wizardModel.setCancelButtonEnabled(new Boolean(newValue));
	}

	/**
	 * Closes the dialog and sets the return code to the integer parameter.
	 * @param code The return code.
	 */
	void close(int code) {
		returnCode = code;
		wizardDialog.dispose();
	}

	private void createButtonPanels() {
		JButton bButton;
		JButton nButton;
		JButton cButton;
		JButton fButton;
		JPanel flowPanel;
		Box box;

		buttonPanel = new JPanel();
		buttonLayout = new CardLayout();
		buttonPanel.setLayout(buttonLayout);

		// Create the buttons with a separator above them, then place them
		// on the east side of the panel with a small amount of space between
		// the back and the next button, and a larger amount of space between
		// the next button and the cancel button.
		
		
		bButton = new JButton(BACK_TEXT, VisualConstants.IMPORT_BACK_ICON);
		nButton = new JButton(NEXT_TEXT, VisualConstants.IMPORT_NEXT_ICON);
		cButton = new JButton(CANCEL_TEXT, VisualConstants.IMPORT_CANCEL_ICON);
		fButton = new JButton(FINISH_TEXT, VisualConstants.IMPORT_FINISH_ICON);

		fButton.setHorizontalTextPosition(SwingConstants.LEFT);
		nButton.setHorizontalTextPosition(SwingConstants.LEFT);
		
		bButton.setEnabled(false);
		nButton.setEnabled(false);
		fButton.setEnabled(false);
		
		bButton.setActionCommand(BACK_BUTTON_ACTION_COMMAND);
		nButton.setActionCommand(NEXT_BUTTON_ACTION_COMMAND);
		cButton.setActionCommand(CANCEL_BUTTON_ACTION_COMMAND);
		fButton.setActionCommand(FINISH_BUTTON_ACTION_COMMAND);

		bButton.addActionListener(wizardController);
		nButton.addActionListener(wizardController);
		cButton.addActionListener(wizardController);
		fButton.addActionListener(wizardController);

		backButtons.put(START_NAVIGATION_PANEL, bButton);
		nextButtons.put(START_NAVIGATION_PANEL, nButton);
		finishButtons.put(START_NAVIGATION_PANEL, fButton);
		cancelButtons.put(START_NAVIGATION_PANEL, cButton);

		flowPanel = new JPanel();
		flowPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
		box = new Box(BoxLayout.X_AXIS);
		flowPanel.add(box);
		buttonPanel.add(flowPanel, START_NAVIGATION_PANEL);

		box.setBorder(new EmptyBorder(new Insets(5, 10, 5, 10)));
		box.add(bButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(nButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(fButton);
		box.add(Box.createHorizontalStrut(30));
		box.add(cButton);

		// Create the navigation panel for back, next and cancel buttons

		bButton = new JButton(BACK_TEXT, VisualConstants.IMPORT_BACK_ICON);
		nButton = new JButton(NEXT_TEXT, VisualConstants.IMPORT_NEXT_ICON);
		cButton = new JButton(CANCEL_TEXT, VisualConstants.IMPORT_CANCEL_ICON);
		
		nButton.setHorizontalTextPosition(SwingConstants.LEFT);

		bButton.setActionCommand(BACK_BUTTON_ACTION_COMMAND);
		nButton.setActionCommand(NEXT_BUTTON_ACTION_COMMAND);
		cButton.setActionCommand(CANCEL_BUTTON_ACTION_COMMAND);

		bButton.addActionListener(wizardController);
		nButton.addActionListener(wizardController);
		cButton.addActionListener(wizardController);

		backButtons.put(NEXT_NAVIGATION_PANEL, bButton);
		nextButtons.put(NEXT_NAVIGATION_PANEL, nButton);
		cancelButtons.put(NEXT_NAVIGATION_PANEL, cButton);

		flowPanel = new JPanel();
		flowPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
		box = new Box(BoxLayout.X_AXIS);
		flowPanel.add(box);
		buttonPanel.add(flowPanel, NEXT_NAVIGATION_PANEL);

		box.setBorder(new EmptyBorder(new Insets(5, 10, 5, 10)));
		box.add(bButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(nButton);
		box.add(Box.createHorizontalStrut(30));
		box.add(cButton);

		// Create the navigation panel for back, finish and cancel buttons

		bButton = new JButton(BACK_TEXT, VisualConstants.IMPORT_BACK_ICON);
		nButton = new JButton(NEXT_TEXT, VisualConstants.IMPORT_NEXT_ICON);
		cButton = new JButton(CANCEL_TEXT, VisualConstants.IMPORT_CANCEL_ICON);
		fButton = new JButton(FINISH_TEXT, VisualConstants.IMPORT_FINISH_ICON);

		fButton.setHorizontalTextPosition(SwingConstants.LEFT);
		nButton.setHorizontalTextPosition(SwingConstants.LEFT);
		
		nButton.setEnabled(false);
		
		bButton.setActionCommand(BACK_BUTTON_ACTION_COMMAND);
		nButton.setActionCommand(NEXT_BUTTON_ACTION_COMMAND);
		cButton.setActionCommand(CANCEL_BUTTON_ACTION_COMMAND);
		fButton.setActionCommand(FINISH_BUTTON_ACTION_COMMAND);

		bButton.addActionListener(wizardController);
		nButton.addActionListener(wizardController);
		cButton.addActionListener(wizardController);
		fButton.addActionListener(wizardController);

		backButtons.put(FINISH_NAVIGATION_PANEL, bButton);
		nextButtons.put(FINISH_NAVIGATION_PANEL, nButton);
		finishButtons.put(FINISH_NAVIGATION_PANEL, fButton);
		cancelButtons.put(FINISH_NAVIGATION_PANEL, cButton);

		flowPanel = new JPanel();
		flowPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
		box = new Box(BoxLayout.X_AXIS);
		flowPanel.add(box);
		buttonPanel.add(flowPanel, FINISH_NAVIGATION_PANEL);

		box.setBorder(new EmptyBorder(new Insets(5, 10, 5, 10)));
		box.add(bButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(nButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(fButton);
		box.add(Box.createHorizontalStrut(30));
		box.add(cButton);

		// Create the navigation panel for back, next, finish and cancel buttons

		bButton = new JButton(BACK_TEXT, VisualConstants.IMPORT_BACK_ICON);
		nButton = new JButton(NEXT_TEXT, VisualConstants.IMPORT_NEXT_ICON);
		cButton = new JButton(CANCEL_TEXT, VisualConstants.IMPORT_CANCEL_ICON);
		fButton = new JButton(FINISH_TEXT, VisualConstants.IMPORT_FINISH_ICON);

		nButton.setHorizontalTextPosition(SwingConstants.LEFT);
		fButton.setHorizontalTextPosition(SwingConstants.LEFT);
		
		bButton.setActionCommand(BACK_BUTTON_ACTION_COMMAND);
		nButton.setActionCommand(NEXT_BUTTON_ACTION_COMMAND);
		cButton.setActionCommand(CANCEL_BUTTON_ACTION_COMMAND);
		fButton.setActionCommand(FINISH_BUTTON_ACTION_COMMAND);

		bButton.addActionListener(wizardController);
		nButton.addActionListener(wizardController);
		cButton.addActionListener(wizardController);
		fButton.addActionListener(wizardController);

		backButtons.put(FINISH_NEXT_NAVIGATION_PANEL, bButton);
		nextButtons.put(FINISH_NEXT_NAVIGATION_PANEL, nButton);
		finishButtons.put(FINISH_NEXT_NAVIGATION_PANEL, fButton);
		cancelButtons.put(FINISH_NEXT_NAVIGATION_PANEL, cButton);

		flowPanel = new JPanel();
		flowPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
		box = new Box(BoxLayout.X_AXIS);
		flowPanel.add(box);
		buttonPanel.add(flowPanel, FINISH_NEXT_NAVIGATION_PANEL);

		box.setBorder(new EmptyBorder(new Insets(5, 10, 5, 10)));
		box.add(bButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(nButton);
		box.add(Box.createHorizontalStrut(10));
		box.add(fButton);
		box.add(Box.createHorizontalStrut(30));
		box.add(cButton);
	}

	/**
	 * This method initializes the components for the wizard dialog: it creates
	 * a JDialog as a CardLayout panel surrounded by a small amount of space on
	 * each side, as well as three buttons at the bottom.
	 */
	private void initComponents() {
		wizardModel.addPropertyChangeListener(this);
		wizardController = new WizardController(this);

		wizardDialog.getContentPane().setLayout(new BorderLayout());
		wizardDialog.addWindowListener(this);

		// Create the outer wizard panel, which is responsible for three
		// buttons:
		//  Next, Back, and Cancel. It is also responsible a JPanel above them that
		//  uses a CardLayout layout manager to display multiple panels in the 
		//  same spot.

		contentPanel = new JPanel();
		contentPanel.setBorder(new EmptyBorder(new Insets(5, 10, 5, 10)));
		
		contentLayout = new CardLayout();
		contentPanel.setLayout(contentLayout);

		JPanel navigationPanel = new JPanel();
		navigationPanel.setLayout(new BorderLayout());

		JSeparator separator = new JSeparator();
		navigationPanel.add(separator, BorderLayout.NORTH);

		createButtonPanels();
		navigationPanel.add(buttonPanel, BorderLayout.EAST);

		wizardDialog.getContentPane().add(navigationPanel,
				java.awt.BorderLayout.SOUTH);
		wizardDialog.getContentPane().add(contentPanel,
				java.awt.BorderLayout.CENTER);

	}

	/**
	 * If the user presses the close box on the dialog's window, treat it
	 * as a cancel.
	 * @param WindowEvent The event passed in from AWT.
	 */

	public void windowClosing(WindowEvent e) {
		returnCode = CANCEL_RETURN_CODE;
	}
}
