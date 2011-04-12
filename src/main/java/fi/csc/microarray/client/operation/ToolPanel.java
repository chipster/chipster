package fi.csc.microarray.client.operation;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.WindowConstants;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;
import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.OperationDefinition.Suitability;
import fi.csc.microarray.client.operation.parameter.ToolParameterPanel;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * The main panel for all tool and parameter choices in
 * the client mainframe.
 * 
 * @author Janne Käki, Petri Klemelä, Aleksi Kallio
 *
 */
public class ToolPanel extends JPanel
							implements ActionListener, PropertyChangeListener {

    // UI texts
	private static final String OPERATION_LIST_TITLE = "Analysis tools";
	private static final String SHOW_PARAMETERS_TEXT = "Show parameters";
	private static final String HIDE_PARAMETERS_TEXT = "Hide parameters";	
	
	// Card layout names
	private static final String TOOLS = "Tools";
	private static final String TOOLS_CATEGORIZED = "Categorized tools";
	private static final String TOOLS_FILTERED = "Filtered tools";
	private static final String PARAMETERS = "Parameters";
	
	// Panel size
	private static final int WHOLE_PANEL_HEIGHT = 240;
	private static final int WHOLE_PANEL_WIDTH = 660;
	
	private JPanel operationPanel;
	private JPanel operationCardPanel;
	private ToolSelectorPanel operationChoicePanel;
	private ToolFilterPanel operationFilterPanel;
	private JTextField searchField;
	private JButton clearSearchButton;
	private JPanel cardPanel;
	private JTextArea detailField = new JTextArea();
	
	private JLabel suitabilityLabel = new JLabel();
	private JButton sourceButton = new JButton("Show tool sourcecode");
	private JButton helpButton = new JButton("More help");
	private JButton parametersButton = new JButton();
	private JButton executeButton = new JButton();
	private JScrollPane detailFieldScroller;
	private boolean isParametersVisible = false;
	
	private OperationDefinition selectedOperationDefinition = null;
	private Operation currentOperation = null;

	
	private ClientApplication application = Session.getSession().getApplication();
	
	/**
	 * Creates a new ToolPanel.
	 * 
	 * @param client The client under whose command this panel is assigned.
	 */
	public ToolPanel(List<ToolCategory> parsedCategories) throws ParseException {
		super(new GridBagLayout());
		this.setPreferredSize(new Dimension(WHOLE_PANEL_WIDTH, WHOLE_PANEL_HEIGHT));
		this.setMinimumSize(new Dimension(0,0));
		
		operationChoicePanel = new ToolSelectorPanel(this, parsedCategories);
		operationFilterPanel = new ToolFilterPanel(this, parsedCategories);
		
		cardPanel = new JPanel(new CardLayout());
		
		detailField.setEditable(false);
		detailField.setLineWrap(true);
		detailField.setWrapStyleWord(true);	
		
		detailFieldScroller = new JScrollPane(detailField);
		detailFieldScroller.setHorizontalScrollBarPolicy(
		        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);		

		sourceButton.setEnabled(false);
        sourceButton.setToolTipText("View Source Code");
		sourceButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {				
				try {
					application.showSourceFor(selectedOperationDefinition.getID());
				} catch (TaskException je) {
					application.reportException(je);
				}
			}			
		});				
				
		helpButton.setEnabled(false);
		helpButton.setToolTipText("More information about this tool");
		helpButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				application.viewHelpFor(selectedOperationDefinition);
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
		executeButton.setEnabled(false);
		
		detailFieldScroller.setBorder(
				BorderFactory.createMatteBorder(1, 0, 0, 0, VisualConstants.TOOL_LIST_BORDER_COLOR));
		
	    // Search bar
        JTabbedPane tabbedPane = new JTabbedPane();        
        tabbedPane.addTab("Microarrays", operationChoicePanel);
        tabbedPane.addTab("NGS", new JPanel());
        tabbedPane.setPreferredSize(new Dimension(400, 0));
        tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		
        JToolBar searchPanel = new JToolBar();
        searchPanel.setLayout(new FlowLayout(FlowLayout.RIGHT, 0, 1));
        // Text field
        searchField = new JTextField(10);
        searchField.addCaretListener(new CaretListener() {
            public void caretUpdate(CaretEvent e) {
                // Show filtered tools
                JTextField field = (JTextField) e.getSource();
                if (field.getText().length() > 0) {
                    field.setBackground(VisualConstants.COLOR_BLUE_LIGHT);
                    if (!clearSearchButton.isAncestorOf(field)) {
                    	field.add(clearSearchButton);
                    }
                	operationFilterPanel.loadFilteredOperations(field.getText());
                    showOperationCard(TOOLS_FILTERED);
                } else {
                	clearSearch();
                }
            }
        });

        // clear button in the search field
        clearSearchButton = new JButton(VisualConstants.CLOSE_FILE_ICON);
        clearSearchButton.setFocusPainted(false);
        clearSearchButton.setContentAreaFilled(false);
        clearSearchButton.setCursor(new Cursor(Cursor.HAND_CURSOR));
        clearSearchButton.setBorder(null);
        clearSearchButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                clearSearch();
            }
        });
        
        // also clear search with esc
        searchField.addKeyListener(new KeyListener() {

			@Override
			public void keyPressed(KeyEvent event) {
				if (event.getKeyCode() == KeyEvent.VK_ESCAPE) {
					clearSearch();
				}
				
			}

			@Override
			public void keyReleased(KeyEvent arg0) {
			}

			@Override
			public void keyTyped(KeyEvent arg0) {
			}
        	
        });
        
        
        
        searchPanel.add(new JLabel(VisualConstants.MAGNIFIER_ICON));
        searchPanel.add(searchField);
        searchField.setLayout(new FlowLayout(FlowLayout.RIGHT, 0, 0));
        searchField.setPreferredSize(new Dimension(100, 22));
        searchPanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1,
                VisualConstants.TOOL_LIST_BORDER_COLOR));
		
		// Operation choice card contains two other cards:
		// operations with categories and filtered operations
		operationPanel = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
	    c.gridx = 0;
        c.gridy = 0;
        c.weightx = 1;
        c.weighty = 0;
        c.gridheight = 1;
        c.insets.set(0,0,0,0);
        c.fill = GridBagConstraints.BOTH;
        searchPanel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
        searchPanel.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);
        operationPanel.add(searchPanel, c);
        operationCardPanel = new JPanel(new CardLayout());
        c.gridy = 1;
        c.weightx = 1;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.insets.set(0,0,0,0);
        operationPanel.add(operationCardPanel, c);
        
        // Tool selection panels inside operation card
        operationCardPanel.add(tabbedPane, TOOLS_CATEGORIZED);
        operationCardPanel.add(operationFilterPanel, TOOLS_FILTERED);
        operationCardPanel.setBorder(BorderFactory.createMatteBorder(1, 0, 0, 0,
                VisualConstants.TOOL_LIST_BORDER_COLOR));

	    // Add operation panel
	    cardPanel.add(operationPanel, TOOLS);
        cardPanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1,
                VisualConstants.TOOL_LIST_BORDER_COLOR));
		
	    // Help and execution panel
		JPanel topLeftPanel = new JPanel(new GridBagLayout());
		c = new GridBagConstraints();
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
		
		// Add everything to the main panel
		c.gridx = 0;
		c.gridy = 0;
		c.gridheight = 3;
		c.weightx = 1;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;
		c.insets.set(0,0,0,0);
		this.add(cardPanel, c);
		c.gridy = 0;
		c.gridx = 1;
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
		System.out.println("FIXME");
//		Session.getSession().getApplication().addClientEventListener(this);		
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
		    executeCurrentOperation();
		}
	}

	private void executeCurrentOperation() {
		// Check if we can run the operation
		Suitability suitability = evaluateSuitability();
		
		if (!suitability.isOk()) {
		    application.showDialog("Check parameters", suitability.toString(), "",
		                           Severity.INFO, true,
		                           DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
		    return;
		}
		
		// Run it	    
		try {
			// we MUST clone the operation, or otherwise results share the same
			// operation as long as it is executed and parameter panel is not closed
			Operation clonedOperation = new Operation(currentOperation);
			application.executeOperation(clonedOperation);
		} catch (MicroarrayException me) {
			throw new RuntimeException(me);
		}
	}

	private void parametersButtonClicked(){
		if(!isParametersVisible){
			showParameterPanel();			
		} else {
			showOperationsPanel();
		}
	}		
	
	/**
	 * Activate a certain card in the left panel.
	 */
	private void showCard(String card) {
	    CardLayout cl = (CardLayout)(cardPanel.getLayout());
        cl.show(cardPanel, card);
	}
	
	/**
     * Activate a certain card in the operation panel.
     * 
     * Currently you can choose between categorized operations and
     * a list of filtered operations.
     */
    private void showOperationCard(String card) {
        CardLayout cl = (CardLayout)(operationCardPanel.getLayout());
        cl.show(operationCardPanel, card);
    }

	private void showParameterPanel() {			
		showParametersTitle(true);
		parametersButton.setText(HIDE_PARAMETERS_TEXT);

		try {
		    // Display parameters for selected operation
			cardPanel.add(new ToolParameterPanel(currentOperation, this),PARAMETERS);
			showCard(PARAMETERS);
		    isParametersVisible = true;

		} catch (MicroarrayException e) {
			application.reportException(e);
		}
	}

	private void showOperationsPanel() {       
        showParametersTitle(false);
        parametersButton.setText(SHOW_PARAMETERS_TEXT);
        isParametersVisible = false;
        showCard(TOOLS);
	}
	
	/**
	 * Sets the title of this panel (shown in the TitledBorder) according to
	 * the received command word (suggesting either a generic "Operations"
	 * title or a more specific one, if an operation is already selected).
	 * The selected dataset will also be taken into account when setting
	 * the title.
	 * 
	 * @param commandWord Either ToolPanel.showOperationsCommand or
	 * 					  ToolPanel.showParametersCommand.
	 */
	private void showParametersTitle(boolean showParametersTitle) {
		String title;
		if (showParametersTitle) {
			// ExecutionItem oper =
			//	operationChoiceView.getSelectedOperation();
			title = OPERATION_LIST_TITLE + 
			  " - " + this.currentOperation.getCategoryName() +
			  " - " + this.currentOperation.getDisplayName();
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
	 * @throws MicroarrayException 
	 */
	public void selectTool(ExecutionItem newSelectedOperationDefinition) {
		
		// no operation selected
		if (newSelectedOperationDefinition == null) {
			clearOperationSelection();
		} 
		
		// operation selected
		else {

			// create the operation
			this.selectedOperationDefinition = (OperationDefinition)newSelectedOperationDefinition;
			try {
				this.currentOperation = new Operation(this.selectedOperationDefinition, application.getSelectionManager().getSelectedDatasAsArray());
			} catch (MicroarrayException e) {
				clearOperationSelection();
				return;
			}

			// show operation info panel
			this.showParametersTitle(false);
			this.showOperationInfoText();
			parametersButton.setEnabled(true);

			// update suitability label and run buttons
			updateSuitability();
		}
	}

	public void updateSuitability() {
		Suitability suitability = evaluateSuitability();
		if (suitability.isOk()) {
			suitabilityLabel.setIcon(VisualConstants.SUITABLE_ICON);
			executeButton.setEnabled(true);
		} else {
			suitabilityLabel.setIcon(VisualConstants.INCOMPATIBLE_ICON);
			executeButton.setEnabled(false);
		}

		suitabilityLabel.setToolTipText(" " + suitability.toString());
	}

	private void clearOperationSelection() {
		this.selectedOperationDefinition = null;
		this.currentOperation = null;

		executeButton.setEnabled(false);
		parametersButton.setEnabled(false);
		suitabilityLabel.setIcon(null);
		suitabilityLabel.setToolTipText("");
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
		if (currentOperation != null) {
			setInfoText(currentOperation.getDescription(), Color.BLACK, true);
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
			detailField.setMargin(new Insets(2,2,2,2+ detailFieldScroller.getVerticalScrollBar().getWidth()));
		} else {
			detailField.setMargin(new Insets(2,2,2,2));
		}
	}
		
	
	/**
	 * Dataset selection changed.
	 */
	public void propertyChange(PropertyChangeEvent dataEvent) {
		if (dataEvent instanceof DatasetChoiceEvent) {
			
            // Reselect operation definition, so suitability is recalculated
			// and parameter values are set to default
            selectTool(selectedOperationDefinition);
            
            // In case parameter panel is open, open operations panel,
            // because we have just nulled the parameter values
            showOperationsPanel();
		}
	}
	
	/**
	 * Reevaluate the suitability of parameter values
	 * and inputs.
	 */
	private Suitability evaluateSuitability() {
		if (currentOperation == null) {
			return Suitability.IMPOSSIBLE;
		}
		
		// Check suitability of parameters and inputs
		Suitability suitability = currentOperation.evaluateSuitabilityFor(
		        application.getSelectionManager().getSelectedDataBeans(), null);
		
		return suitability;
	}

	public void runSelectedOperation() {
		if (evaluateSuitability().isOk()) {
			this.executeCurrentOperation();
		}
	}

	private void clearSearch() {
		searchField.setText("");
		searchField.setBackground(Color.WHITE);
		searchField.remove(clearSearchButton);
		operationChoicePanel.deselectTool();
		showOperationCard(TOOLS_CATEGORIZED);
	}
	
	public static void main(String[] args) throws Exception {
		ToolPanel toolPanel = new ToolPanel(new LinkedList<ToolCategory>());
		JFrame frame = new JFrame();
		frame.add(toolPanel);
		frame.setSize(1000, 400);
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}
}
