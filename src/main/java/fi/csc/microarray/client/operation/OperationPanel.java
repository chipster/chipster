package fi.csc.microarray.client.operation;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Collection;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;

import org.apache.log4j.Logger;

import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition.Suitability;
import fi.csc.microarray.client.operation.parameter.ToolParameterPanel;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * The main panel for all operation, parameter and visualization choices in
 * the client mainframe.
 * 
 * @author Janne KÃ¤ki, Petri KlemelÃ¤
 *
 */
public class OperationPanel extends JPanel
							implements ActionListener, PropertyChangeListener {

	private static final String OPERATION_LIST_TITLE = "Analysis tools";
	private static final String SHOW_PARAMETERS_TEXT = "Show parameters";
	private static final String HIDE_PARAMETERS_TEXT = "Hide parameters";	
	
	private static final String OPERATIONS = "Operations";
	private static final String PARAMETERS = "Parameters";

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(OperationPanel.class);

	private static final int WHOLE_PANEL_HEIGHT = 240;
	private static final int WHOLE_PANEL_WIDTH= 660;
	
	private OperationChoicePanel operationChoicePanel = null;
	private JPanel cardPanel = new JPanel();
	private JTextArea detailField = new JTextArea();
	
	private JLabel suitabilityLabel = new JLabel();
	private JButton sourceButton = new JButton("Show tool sourcecode");
	private JButton helpButton = new JButton("More help");
	private JButton parametersButton = new JButton();
	private JButton executeButton = new JButton();
	private JScrollPane detailFieldScroller;
	private boolean isParametersVisible = false;
	
	private ExecutionItem chosenOperation = null;
	private ClientApplication application = Session.getSession().getApplication();

	
	/**
	 * Creates a new OperationPanel.
	 * 
	 * @param client The client under whose command this panel is assigned.
	 */
	public OperationPanel(Collection<OperationCategory> parsedCategories) throws ParseException {
		super(new GridBagLayout());
		this.setPreferredSize(new Dimension(WHOLE_PANEL_WIDTH, WHOLE_PANEL_HEIGHT));
		this.setMinimumSize(new Dimension(0,0));
		
		operationChoicePanel = new OperationChoicePanel(this, parsedCategories);		
		
		cardPanel = new JPanel(new CardLayout());
		
		detailField.setEditable(false);
		detailField.setLineWrap(true);
		detailField.setWrapStyleWord(true);	
		
		detailFieldScroller = new JScrollPane(detailField);
		detailFieldScroller.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);		

		sourceButton.setEnabled(false);
        sourceButton.setToolTipText("View Source Code");
		sourceButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {				
				try {
					application.showSourceFor(SADLParser.generateOperationIdentifier(chosenOperation.getCategoryName(), chosenOperation.getName()));
				} catch (TaskException je) {
					application.reportException(je);
				}
			}			
		});				
				
		helpButton.setEnabled(false);
		helpButton.setToolTipText("More information about this tool");
		helpButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				OperationDefinition od = chosenOperation instanceof OperationDefinition ? (OperationDefinition) chosenOperation : ((Operation)chosenOperation).getDefinition();
				application.viewHelpFor(od);
			}			
		});
				
		parametersButton.addActionListener(this);
		parametersButton.setEnabled(false);
		parametersButton.setText(SHOW_PARAMETERS_TEXT);				

		suitabilityLabel.setPreferredSize(new Dimension(
				VisualConstants.INCOMPATIBLE_ICON.getIconHeight(),
				VisualConstants.INCOMPATIBLE_ICON.getIconHeight()));		
		
		executeButton.setIcon(VisualConstants.DOUBLE_FORWARD_ICON);
		executeButton.setDisabledIcon(VisualConstants.DOUBLE_FORWARD_BW_ICON);
		executeButton.setText("<html><b>Run</b></html>");
		executeButton.setHorizontalAlignment(SwingConstants.CENTER);
		executeButton.setHorizontalTextPosition(SwingConstants.LEFT);
		executeButton.setToolTipText("Run selected operation for selected datasets");
		executeButton.addActionListener(this);
		executeButton.setName("executeButton");
		setExecuteButtonEnabled(false);
		
		detailFieldScroller.setBorder(
				BorderFactory.createMatteBorder(1, 0, 0, 0, VisualConstants.OPERATION_LIST_BORDER_COLOR));
		
		cardPanel.add(operationChoicePanel,OPERATIONS);
		
		JPanel topLeftPanel = new JPanel(new GridBagLayout());				
		
		GridBagConstraints c = new GridBagConstraints();
		
		c.gridx = 0;
		c.gridy = 0;
		c.gridheight = 1;
		c.weightx = 0;
		c.weighty = 0;				
		c.fill = GridBagConstraints.NONE;
		c.insets.set(0,10,0,10);
		
		topLeftPanel.add(suitabilityLabel, c);
		c.gridx++;
		c.weightx = 1;
		c.weighty = 1;				
		c.fill = GridBagConstraints.BOTH;
		c.insets.set(0,0,0,0);
		topLeftPanel.add(parametersButton,c);		
		
		JPanel topPanel = new JPanel(new GridLayout(1,2));
		topPanel.add(topLeftPanel);
		topPanel.add(executeButton);
		JPanel bottomPanel = new JPanel(new GridLayout(1,2));
		bottomPanel.add(helpButton);
		bottomPanel.add(sourceButton);
		
		
		c.gridx = 0;
		c.gridy = 0;
		c.gridheight = 3;
		c.weightx = 1;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.insets.set(0,0,0,0);
		this.add(cardPanel, c);
		c.gridx++;
		c.gridwidth = 1;
		c.gridheight = 1;
		c.weightx = 0;
		c.weighty = 0;
		this.add(topPanel,c);		
		c.gridx = 1;
		c.gridy++;		
		c.gridheight = 1;
		c.weightx = 0;
		c.weighty = 1;
		this.add(detailFieldScroller, c);
		c.gridy++;		
		c.weightx = 0;
		c.weighty = 0;
		this.add(bottomPanel, c);		

		// start listening
		Session.getSession().getApplication().addPropertyChangeListener(this);
		
	}
	
	public Vector<Component> getFocusComponents(){
				
		Vector<Component> order = new Vector<Component>();
		order.addAll(operationChoicePanel.getFocusComponents());
		//No easy way to add parametercomponents
		//order.add(parametersButton);
		order.add(executeButton);				
		return order;
	}
	
	/**
	 * The method for listening the top right corner arrow button.
	 * 
	 * @param e The action event.
	 */
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == parametersButton) {
			parametersButtonClicked();
		} else if (e.getSource() == executeButton ) {
			if (chosenOperation instanceof OperationDefinition) {
				application.executeOperation((OperationDefinition)chosenOperation, null);
			} else {				
				try {
					// we MUST clone the operation, or otherwise results share the same
					// operation as long as it is executed and parameter panel is not closed
					Operation clonedOperation = new Operation((Operation)chosenOperation);
					application.executeOperation(clonedOperation);
				} catch (MicroarrayException me) {
					throw new RuntimeException(me);
				}
			}
		}
	}

	private void parametersButtonClicked(){
		if(!isParametersVisible){
			showParameterPanel();			
		} else {
			showOperationsPanel();
		}
	}		
	
	private void setExecuteButtonEnabled(boolean enabled){
		if (enabled){
			executeButton.setEnabled(true);
		} else {
			executeButton.setEnabled(false);
		}
	}
	

	private void showParameterPanel() {			
			showParametersTitle(true);
			parametersButton.setText(HIDE_PARAMETERS_TEXT);

			try {

				if (chosenOperation instanceof Operation) {
					chosenOperation = (Operation)chosenOperation;

				} else if (chosenOperation instanceof OperationDefinition) {
					chosenOperation = new Operation((OperationDefinition)chosenOperation, application.getSelectionManager().getSelectedDatasAsArray());

				} else {
					throw new RuntimeException("wrong type: " + chosenOperation.getClass().getSimpleName());				
				}
				
				cardPanel.add(new ToolParameterPanel((Operation)chosenOperation, this),PARAMETERS);
				CardLayout cl = (CardLayout)(cardPanel.getLayout());
			    cl.show(cardPanel, PARAMETERS);
			    isParametersVisible = true;

			} catch (MicroarrayException e) {
				application.reportException(e);
			}		
	}

	private void showOperationsPanel() {
		showParametersTitle(false);
		parametersButton.setText(SHOW_PARAMETERS_TEXT);
		CardLayout cl = (CardLayout)(cardPanel.getLayout());
		 cl.show(cardPanel, OPERATIONS);
		 isParametersVisible = false;
	}
	
	/**
	 * Sets the title of this panel (shown in the TitledBorder) according to
	 * the received command word (suggesting either a generic "Operations"
	 * title or a more specific one, if an operation is already selected).
	 * The selected dataset will also be taken into accound when setting
	 * the title.
	 * 
	 * @param commandWord Either OperationPanel.showOperationsCommand or
	 * 					  OperationPanel.showParametersCommand.
	 */
	private void showParametersTitle(boolean showParametersTitle) {
		String title;
		if (showParametersTitle) {
			// ExecutionItem oper =
			//	operationChoiceView.getSelectedOperation();
			title = OPERATION_LIST_TITLE + 
			  " - " +this.chosenOperation.getCategoryName() +
			  " [" + this.chosenOperation.getName() + "]";
		} else {
			title = OPERATION_LIST_TITLE;
		}
		((SimpleInternalFrame)this.getParent()).setTitle(title);
		this.repaint();
	}
		
	/**
	 * Shows the "details" (that is, description text and suitability
	 * evaluation) for the given operation. Suitability is evaluated for
	 * the currently chosen dataset.
	 * 
	 * @param operation The operation (or operation definition, or workflow)
	 * 					whose details are to be shown.
	 */
	public void selectOperation(ExecutionItem operation) {
		
		// update source button and operation description text
		
		this.chosenOperation = operation;
		this.showParametersTitle(false);
		
		this.showOperationInfoText();
		
		// update suitability label and action buttons
		if (operation != null) {
			parametersButton.setEnabled(true);
			
			Suitability suitability = evaluateSuitability();
			if (suitability.isImpossible()) {
				makeButtonsEnabled(false);
			} else {
				makeButtonsEnabled(true);
			}
			
			if( suitability.isImpossible()){
				suitabilityLabel.setIcon(VisualConstants.INCOMPATIBLE_ICON);
			} else if( suitability.isOk()){
				suitabilityLabel.setIcon(VisualConstants.SUITABLE_ICON);
			} else if( suitability.isImpossible()){
				suitabilityLabel.setIcon(VisualConstants.SUITABILITY_WARNING_ICON);
			}
					
			suitabilityLabel.setToolTipText(" " + suitability.toString());
		} else {
			makeButtonsEnabled(false);
			parametersButton.setEnabled(false);
			suitabilityLabel.setIcon(null);
			suitabilityLabel.setToolTipText("");
		}
	}
	
	/**
	 * Sets new info text and updates text margins if the vertical scrollbar 
	 * appers
	 * 
	 * @param text
	 * @param color
	 * @param enable
	 */
	public void setInfoText(String text, Color color, boolean enable) {
		detailField.setForeground(color);
		detailField.setText(text);
		
		// Increases text margins if the scrollbar appears
		updateDetailFieldMargin();
	}
	
	public void showOperationInfoText() {
		if (chosenOperation != null) {
			setInfoText(chosenOperation.getDescription(), Color.BLACK, true);
			sourceButton.setEnabled(true);
			helpButton.setEnabled(true);
		} else {
			setInfoText("", Color.BLACK, false);
			sourceButton.setEnabled(false);
			helpButton.setEnabled(false);
		}
	}

	/**
	 * Increases the text margin if the scrollbar appers. The problem is that 
	 * word wrapping is done before scrollbar appers and some of the text is left 
	 * behind scrollbar.
	 * 
	 * @author mkoski
	 */
	private void updateDetailFieldMargin(){
		// TODO I think the word wrapping problem should be fixed in a some better way
		// than increasing the margin
		
		// Gets the scroller and sets more margin if the scrollbar is visible
		if(detailFieldScroller.getVerticalScrollBar().isVisible()){
			logger.debug("vertical scrollbar is visible");
			detailField.setMargin(new Insets(2,2,2,2+ detailFieldScroller.getVerticalScrollBar().getWidth()));
		} else {
			logger.debug("vertical scrollbar is not visible");
			detailField.setMargin(new Insets(2,2,2,2));
		}
	}
		
	/**
	 * Enables (or disables) action - the two buttons in the top right corner,
	 * to be exact. This prevents from enabling an operation without selecting
	 * a dataset first.
	 * 
	 * @param enabled Whether action is to be enabled or not.
	 */
	public void enableAction(boolean enabled) {
		if (!evaluateSuitability().isImpossible()) {
			makeButtonsEnabled(enabled);			
		}
		if(!enabled){
			suitabilityLabel.setIcon(null);
		}
	}
	
	private void makeButtonsEnabled(boolean enabled) {
		if(!enabled){
			this.showOperationsPanel();
		}
		//parametersButton.setEnabled(enabled);
		setExecuteButtonEnabled(enabled);
	}
	
	public void propertyChange(PropertyChangeEvent dataEvent) {
		if(dataEvent instanceof DatasetChoiceEvent) {
			logger.debug("chosen data " +  application.getSelectionManager().getSelectedDataBean() + " (possible one among many)");
			if (application.getSelectionManager().getSelectedDataBean() != null) {
				showOperationsPanel();
			}
		}
	}
	
	private Suitability evaluateSuitability() {
		if (chosenOperation == null) {
			return Suitability.IMPOSSIBLE;
		}
		
		return chosenOperation.evaluateSuitabilityFor(application.getSelectionManager().getSelectedDataBeans());
	}
}
