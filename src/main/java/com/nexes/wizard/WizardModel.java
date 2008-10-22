package com.nexes.wizard;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.HashMap;


/**
 * The model for the Wizard component, which tracks the text and enabled state
 * of each of the buttons, as well as the current panel that is displayed. Note that 
 * the model, in its current form, is not intended to be subclassed. 
 */


public class WizardModel {

	/**
	 * Identification string for the current panel.
	 */
	public static final String CURRENT_PANEL_DESCRIPTOR_PROPERTY = "currentPanelDescriptorProperty";

	/**
	 * Property identification String for the Back button's text
	 */
	public static final String BACK_BUTTON_TEXT_PROPERTY = "backButtonTextProperty";

	/**
	 * Property identification String for the Back button's enabled state
	 */
	public static final String BACK_BUTTON_ENABLED_PROPERTY = "backButtonEnabledProperty";

	/**
	 * Property identification String for the Next button's text
	 */
	public static final String NEXT_BUTTON_TEXT_PROPERTY = "nextButtonTextProperty";

	/**
	 * Property identification String for the Next button's enabled state
	 */
	public static final String NEXT_BUTTON_ENABLED_PROPERTY = "nextButtonEnabledProperty";

	/**
	 * Property identification String for the Finish button's text
	 */
	public static final String FINISH_BUTTON_TEXT_PROPERTY = "finishButtonTextProperty";

	/**
	 * Property identification String for the Finish button's enabled state
	 */
	public static final String FINISH_BUTTON_ENABLED_PROPERTY = "finishButtonEnabledProperty";

	/**
	 * Property identification String for the Cancel button's text
	 */
	public static final String CANCEL_BUTTON_TEXT_PROPERTY = "cancelButtonTextProperty";

	/**
	 * Property identification String for the Cancel button's enabled state
	 */
	public static final String CANCEL_BUTTON_ENABLED_PROPERTY = "cancelButtonEnabledProperty";

	private WizardPanelDescriptor currentPanel;

	private HashMap<Object, WizardPanelDescriptor> panelHashmap;

	Object backButtonText;

	Object nextButtonText;

	Object finishButtonText;

	Object cancelButtonText;

	private Boolean backButtonEnabled;

	private Boolean nextButtonEnabled;

	private Boolean finishButtonEnabled;

	private Boolean cancelButtonEnabled;

	private PropertyChangeSupport propertyChangeSupport;

	/**
	 * Default constructor.
	 */
	public WizardModel() {

		panelHashmap = new HashMap<Object, WizardPanelDescriptor>();

		propertyChangeSupport = new PropertyChangeSupport(this);
	}

	/**
	 * Returns the currently displayed WizardPanelDescriptor.
	 * @return The currently displayed WizardPanelDescriptor
	 */
	WizardPanelDescriptor getCurrentPanelDescriptor() {
		return currentPanel;
	}

	/**
	 * Registers the WizardPanelDescriptor in the model using the Object-identifier specified.
	 * @param id Object-based identifier
	 * @param descriptor WizardPanelDescriptor that describes the panel
	 */	
	void registerPanel(Object id, WizardPanelDescriptor descriptor) {

		//  Place a reference to it in a hashtable so we can access it later
		//  when it is about to be displayed.

		panelHashmap.put(id, descriptor);

	}

	/**
	 * Sets the current panel to that identified by the Object passed in.
	 * @param id Object-based panel identifier
	 * @return boolean indicating success or failure
	 */
	boolean setCurrentPanel(Object id) {

		//  First, get the hashtable reference to the panel that should
		//  be displayed.

		WizardPanelDescriptor nextPanel = (WizardPanelDescriptor) panelHashmap
				.get(id);

		//  If we couldn't find the panel that should be displayed, return
		//  false.

		if (nextPanel == null)
			throw new WizardPanelNotFoundException();

		WizardPanelDescriptor oldPanel = currentPanel;
		currentPanel = nextPanel;

		if (oldPanel != currentPanel)
			firePropertyChange(CURRENT_PANEL_DESCRIPTOR_PROPERTY, oldPanel,
					currentPanel);

		return true;

	}

	Object getBackButtonText() {
		return backButtonText;
	}

	void setBackButtonText(Object newText) {

		Object oldText = backButtonText;
		if (!newText.equals(oldText)) {
			backButtonText = newText;
			firePropertyChange(BACK_BUTTON_TEXT_PROPERTY, oldText, newText);
		}
	}

	Object getNextButtonText() {
		return nextButtonText;
	}

	void setNextButtonText(Object newText) {

		Object oldText = nextButtonText;
		if (!newText.equals(oldText)) {
			nextButtonText = newText;
			firePropertyChange(NEXT_BUTTON_TEXT_PROPERTY, oldText, newText);
		}
	}

	Object getFinishButtonText() {
		return finishButtonText;
	}

	void setFinishButtonText(Object newText) {

		Object oldText = finishButtonText;
		if (!newText.equals(oldText)) {
			finishButtonText = newText;
			firePropertyChange(FINISH_BUTTON_TEXT_PROPERTY, oldText, newText);
		}
	}

	Object getCancelButtonText() {
		return cancelButtonText;
	}

	void setCancelButtonText(Object newText) {

		Object oldText = cancelButtonText;
		if (!newText.equals(oldText)) {
			cancelButtonText = newText;
			firePropertyChange(CANCEL_BUTTON_TEXT_PROPERTY, oldText, newText);
		}
	}

	Boolean getBackButtonEnabled() {
		return backButtonEnabled;
	}

	void setBackButtonEnabled(Boolean newValue) {

		Boolean oldValue = backButtonEnabled;
		if (newValue != oldValue) {
			backButtonEnabled = newValue;
			firePropertyChange(BACK_BUTTON_ENABLED_PROPERTY, oldValue, newValue);
		}
	}

	Boolean getNextButtonEnabled() {
		return nextButtonEnabled;
	}

	void setNextButtonEnabled(Boolean newValue) {

		Boolean oldValue = nextButtonEnabled;
		if (newValue != oldValue) {
			nextButtonEnabled = newValue;
			firePropertyChange(NEXT_BUTTON_ENABLED_PROPERTY, oldValue, newValue);
		}
	}

	Boolean getFinishButtonEnabled() {
		return finishButtonEnabled;
	}

	void setFinishButtonEnabled(Boolean newValue) {

		Boolean oldValue = finishButtonEnabled;
		if (newValue != oldValue) {
			finishButtonEnabled = newValue;
			firePropertyChange(FINISH_BUTTON_ENABLED_PROPERTY, oldValue,
					newValue);
		}
	}

	Boolean getCancelButtonEnabled() {
		return cancelButtonEnabled;
	}

	void setCancelButtonEnabled(Boolean newValue) {

		Boolean oldValue = cancelButtonEnabled;
		if (newValue != oldValue) {
			cancelButtonEnabled = newValue;
			firePropertyChange(CANCEL_BUTTON_ENABLED_PROPERTY, oldValue,
					newValue);
		}
	}

	public void addPropertyChangeListener(PropertyChangeListener p) {
		propertyChangeSupport.addPropertyChangeListener(p);
	}

	public void removePropertyChangeListener(PropertyChangeListener p) {
		propertyChangeSupport.removePropertyChangeListener(p);
	}

	protected void firePropertyChange(String propertyName, Object oldValue,
			Object newValue) {
		propertyChangeSupport.firePropertyChange(propertyName, oldValue,
				newValue);
	}

}
