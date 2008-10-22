package fi.csc.microarray.client.dataimport.table;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;
import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.ToolBarComponentFactory;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.dataimport.ColumnType;
import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.events.ChipCountChangeEvent;
import fi.csc.microarray.client.dataimport.events.ChipNumberChangedEvent;
import fi.csc.microarray.client.dataimport.events.ColumnTitlesChangedEvent;
import fi.csc.microarray.client.dataimport.events.ColumnTypeChangeListener;
import fi.csc.microarray.client.dataimport.events.ColumnTypeChangedEvent;
import fi.csc.microarray.client.dataimport.events.ConversionModelChangeListener;
import fi.csc.microarray.client.dataimport.events.DecimalSeparatorChangedEvent;
import fi.csc.microarray.client.dataimport.events.DelimiterChangedEvent;
import fi.csc.microarray.client.dataimport.events.FooterChangedEvent;
import fi.csc.microarray.client.dataimport.events.HeaderChangedEvent;
import fi.csc.microarray.client.dataimport.events.InputFileChangedEvent;
import fi.csc.microarray.client.dataimport.events.TitleRowChangedEvent;

/**
 * Class for import table (the right side of the main split pane)
 * 
 * @author mkoski klemela
 *
 */
public class TableInternalFrame extends SimpleInternalFrame implements
		ActionListener, ChangeListener, ConversionModelChangeListener, 
		ColumnTypeChangeListener{

	private static final String NOT_ENOUGH_CHIPS_TOOLTIP = "Number of selections of this column type is bigger than chip count";
	private static final String CHIP_NUMBER_NOT_SET_TOOLTIP = "Select the chip number";
	
	private static final Dimension SPINNER_SIZE = new Dimension(90,20);
	
	private JToolBar toolBarFirstStep;
	private JToolBar toolBarSecondStep;
	private ImportPreviewTable table;
	private ImportScreen screen;
	private JScrollPane tableScroller;
	private JButton toTopButton;
	private JButton toBottomButton;
	private JButton resetButtonFirstStep;
	private JToggleButton markHeaderButton;
	private JToggleButton markFooterButton;
	private JToggleButton markTitleButton;	
	private List<JToggleButton> markHeaderFooterButtons;
	private JSpinner headerSpinner;
	private JSpinner footerSpinner;
	private JLabel showingColumnsLabel;

	/**
	 * Toolbar buttons for marking sample, flag etc columns.
	 * This map connects button to the column type that the 
	 * button represents.
	 */
	private Map<JToggleButton,ColumnType> markColumnsButtons;
	private JButton resetButtonSecondStep;

	public TableInternalFrame(ImportScreen screen) {
		super("Import data");
		this.screen = screen;
	}

	/**
	 * Initializes the first step
	 *
	 */
	public void initializeFirstStep() {
		this.updateTitle();
		
		// Set limits to the first step (unlimited rows and 4 data column (the first 
		// column is for row numbers, so that's why it is one bigger))
		screen.getConversionModel().setLimits(Integer.MAX_VALUE, 6);
				
		this.setToolBar(getToolbarFirstStep());
		this.setContent(getContentPanel());
		this.updateHilighter();
		
	}
	
	private void updateTitle() {
		if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
			this.setTitle("Select rows (" + screen.getConversionModel().getInputFileName() + ")");
		} else {
			this.setTitle("Select columns (" + screen.getConversionModel().getInputFileName() + ")");
		}
	}

	private JPanel getContentPanel() {
		JPanel contentPanel = new JPanel(new BorderLayout());
		contentPanel.add(getShowingColumnsLabel(), BorderLayout.NORTH);
		contentPanel.add(getTableScroller(), BorderLayout.CENTER);
		return contentPanel;
	}

	private JLabel getShowingColumnsLabel() {
		if(showingColumnsLabel == null){
			showingColumnsLabel = new JLabel();
			updateShowingColumnsLabel();
			return showingColumnsLabel;
		} else {
			return showingColumnsLabel;
		}
	}

	/**
	 * Sets the text for "Showing columns/rows ... of ..." label
	 *
	 */
	public void updateShowingColumnsLabel(){
		if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
			// First step
			
			// -1 because the column count includes the row number column
			int cols = screen.getConversionModel().getLimitedColumnCount() - 1;
			int colsTotal = screen.getConversionModel().getUnlimitedColumnCount() - 1;
			
			// If data is not set yet, show zero (0) rather than minus one (-1)
			if(cols == -1 && colsTotal == -1){
				cols = 0;
				colsTotal = 0;
			}
			
			showingColumnsLabel.setText("Showing columns " + cols + " of " + colsTotal);
		} else {
			// Second step
			int rows = screen.getConversionModel().getLimitedRowCount();
			int rowsTotal = screen.getConversionModel().getUnlimitedRowCount();
			showingColumnsLabel.setText("Showing rows " + rows + " of " + rowsTotal);
		}
	}
	
	/**
	 * Initializes the second step
	 *
	 */
	public void initializeSecondStep() {
		this.updateTitle();
		
		// Sets limits to the second step (100 rows, unlimited columns)
		screen.getConversionModel().setLimits(100, Integer.MAX_VALUE);
		
		this.setToolBar(getToolbarSecondStep());
		this.setContent(getContentPanel());
		
		this.updateHilighter();
	}	

	private void updateHilighter() {
		boolean rowHilight = false;
		boolean columnHilight = false;
		
		if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
			rowHilight = 
				markFooterButton.isSelected() ||
				markHeaderButton.isSelected() ||
				markTitleButton.isSelected();
			
		} else if(screen.getCurrentStep() == ImportScreen.Step.SECOND){
			columnHilight = true;
		}
		
		table.setColumnHighlight(columnHilight);
		table.setRowHighlight(rowHilight);
		
	}

	/**
	 * Gets toolbar for the first step
	 * 
	 * @return
	 */
	private JToolBar getToolbarFirstStep() {
		if (toolBarFirstStep == null) {
			toolBarFirstStep = new JToolBar();
			//Original layout manager streches the spinners
			toolBarFirstStep.setLayout(new GridBagLayout());
			toolBarFirstStep.putClientProperty(Options.HEADER_STYLE_KEY,
					HeaderStyle.SINGLE);
			
			markHeaderButton = ToolBarComponentFactory.createToggleButton(
					"Mark header", VisualConstants.IMPORT_HEADER_ICON, true, true);
			
			markHeaderButton.addActionListener(this);
			markFooterButton = ToolBarComponentFactory.createToggleButton(
					"Mark footer", VisualConstants.IMPORT_FOOTER_ICON, true, true);
		
			markFooterButton.addActionListener(this);
			markTitleButton = ToolBarComponentFactory.createToggleButton(
					"Mark title row", VisualConstants.IMPORT_TITLE_ICON, true,true);
			
			markTitleButton.addActionListener(this);
			headerSpinner = ToolBarComponentFactory.createSpinner();
			headerSpinner.addChangeListener(this);			 
			// The second parameter of NumberEditor constructor is number format
			// The string "0" means that it is a digit without any thousand separators
			// or decimals
			headerSpinner.setEditor(new JSpinner.NumberEditor(headerSpinner, "0"));
			headerSpinner.setValue(0);
			headerSpinner.setPreferredSize(SPINNER_SIZE);
			
			footerSpinner = ToolBarComponentFactory.createSpinner();
			footerSpinner.addChangeListener(this);
			footerSpinner.setEditor(new JSpinner.NumberEditor(footerSpinner, "0"));	
			footerSpinner.setValue(Integer.MAX_VALUE);
			footerSpinner.setPreferredSize(SPINNER_SIZE);
			
			resetButtonFirstStep = ToolBarComponentFactory.createButton("Reset", VisualConstants.IMPORT_RESET_ICON, true, false);
			resetButtonFirstStep.addActionListener(this);

			GridBagConstraints c = new GridBagConstraints();
			c.weightx = 0;
			c.anchor = GridBagConstraints.WEST;
			c.fill = GridBagConstraints.VERTICAL;
			toolBarFirstStep.add(markHeaderButton,c);
			toolBarFirstStep.add(headerSpinner,c);
			toolBarFirstStep.addSeparator();
			toolBarFirstStep.add(markFooterButton,c);
			toolBarFirstStep.add(footerSpinner,c);
			toolBarFirstStep.addSeparator();
			toolBarFirstStep.add(markTitleButton,c);
			c.weightx = 1;
			c.fill = GridBagConstraints.BOTH;
			toolBarFirstStep.add(new JLabel(),c);			
			c.weightx = 0;
			c.fill = GridBagConstraints.VERTICAL;
			toolBarFirstStep.add(resetButtonFirstStep, c);

			// Set buttons to list
			markHeaderFooterButtons = new ArrayList<JToggleButton>();
			markHeaderFooterButtons.add(markHeaderButton);
			markHeaderFooterButtons.add(markFooterButton);
			markHeaderFooterButtons.add(markTitleButton);
		}
		return toolBarFirstStep;
	}

	/**
	 * Gets toolbar for the second step
	 * @return
	 */
	private JToolBar getToolbarSecondStep() {
		if (toolBarSecondStep == null){
			toolBarSecondStep = new JToolBar();
			toolBarSecondStep.setLayout(new GridLayout(1,ColumnType.values().length));
			toolBarSecondStep.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);
			
			markColumnsButtons = new HashMap<JToggleButton,ColumnType>();
			for(int i = 0; i < ColumnType.values().length; i++){
				if(ColumnType.values()[i].equals(ColumnType.ROW_NUMBER)){
					// skip row number
					continue;
				}
				JToggleButton button;

				button = ToolBarComponentFactory.createToggleButton(
						ColumnType.values()[i].toString(), //text
						null,								//icon
						false,		//left border
						true);		//right border
				
				
				if(i == 0){
					button.setSelected(true);
				}
				
				button.addActionListener(this);		
				
				markColumnsButtons.put(button,ColumnType.values()[i]);
				toolBarSecondStep.add(button);
			}			
			
			resetButtonSecondStep = ToolBarComponentFactory.createButton("Reset", VisualConstants.IMPORT_RESET_ICON, true, true);
			resetButtonSecondStep.addActionListener(this);
			toolBarSecondStep.add(resetButtonSecondStep);
			
			this.setToolBar(toolBarSecondStep);
		}
		return toolBarSecondStep;
	}

	/**
	 * Selects a button from the list and deselects other 
	 * buttons in the list
	 * @param selected the selected button
	 * @param buttons button list
	 */
	private static void selectButton(JToggleButton selected,
			Collection<JToggleButton> buttons) {
		for (JToggleButton button : buttons) {
			if (button == selected) {
				button.setSelected(true);
			} else {
				button.setSelected(false);
			}
		}
	}

	/**
	 * Gets table scroller
	 * 
	 * @return
	 */
	public JScrollPane getTableScroller() {
				
		tableScroller = new JScrollPane(getTable());
		if(this.screen.getCurrentStep() == ImportScreen.Step.FIRST){
			addCornerComponents();
		}
		return tableScroller;
	}
	
	public void addCornerComponents(){

		toTopButton = new JButton(VisualConstants.TO_TOP_ICON); 
		toBottomButton = new JButton(VisualConstants.TO_BOTTOM_ICON);

		toTopButton.setToolTipText("Go to beginning of the table");
		toBottomButton.setToolTipText("Go to end of the table");

		toTopButton.addActionListener(this);
		toBottomButton.addActionListener(this);

		toTopButton.setBorder(BorderFactory.createMatteBorder(0,1,1,0, 
				UIManager.getColor("ScrollBar.darkShadow")));
		toBottomButton.setBorder(BorderFactory.createMatteBorder(0,1,0,0, 
				UIManager.getColor("ScrollBar.darkShadow")));

		tableScroller.setCorner(JScrollPane.UPPER_RIGHT_CORNER, toTopButton);			
		tableScroller.setCorner(JScrollPane.LOWER_RIGHT_CORNER, toBottomButton);

		tableScroller.setHorizontalScrollBarPolicy(
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		//No scrolling needed, just to make space for the to bottom button
		tableScroller.getHorizontalScrollBar().setEnabled(false);
	}

	/**
	 * Gets the table
	 * 
	 * @return
	 */
	public ImportPreviewTable getTable() {
		if (table == null) {
			table = new ImportPreviewTable(screen, this);
			table.setCellSelectionEnabled(false);
			screen.getConversionModel().addConversionChangeListener(table);
			return table;
		} else {
			return table;
		}
	}
	
	/**
	 * Updates all chip count comboboxes. This should be done when the 
	 * chip count changes for example.
	 *
	 */
	public void updateAllChipCountComboBoxes(){
		int columnCount = screen.getColumnTypeManager().getColumnCount();

		// Iterater through the columns (first column (0) is null)
		for(int columnIndex = 1; columnIndex < columnCount; columnIndex++){
			updateChipCountComboBox(columnIndex);
		}
	}
	
	/**
	 * Updates chip count combobox.
	 * 
	 * @param columnIndex
	 */
	public void updateChipCountComboBox(int columnIndex){
		int chipCount = screen.getColumnTypeManager().getChipCount();
		
		Object[] items;
		Object selectedItem = null;

		// Keep the list empty if the type of column is UNUSED
		if(screen.getColumnTypeManager().getColumnType(columnIndex).equals(ColumnType.UNUSED_LABEL)
				|| screen.getColumnTypeManager().getColumnType(columnIndex).equals(ColumnType.ANNOTATION_LABEL)
				|| screen.getColumnTypeManager().getColumnType(columnIndex).equals(ColumnType.IDENTIFIER_LABEL)){
			items = new Integer[0];
		} else {
			items = new Integer[chipCount];
			for(int chipNum = 1; chipNum <= chipCount; chipNum++){
				items[chipNum-1] = chipNum;
				if(chipNum == screen.getColumnTypeManager().getColumnChipNumber(columnIndex)){
					selectedItem = chipNum;
				}
			}
		}

		ComboBoxModel model = new DefaultComboBoxModel(items);
		model.setSelectedItem(selectedItem);
				
		table.getHeaderRenderer(columnIndex).getCombo().setModel(model);
		table.getHeaderRenderer(columnIndex).getCombo().setEnabled(model.getSize() > 0);
		table.getHeaderRenderer(columnIndex).update();
	}
	
	/**
	 * Updates all column title labels. This should be done for example when the 
	 * chip count changed.
	 *
	 */
	public void updateAllColumnTitleLabels(){
		// Do not update label 0, because it is null
		for(int i = 1; i < table.getColumnCount(); i++){
			updateColumnTitleLabel(i);
		}
	}
	
	/**
	 * Updates column title label
	 * @param columnIndex
	 */
	public void updateColumnTitleLabel(int columnIndex) {
		if(screen.getColumnTypeManager().getColumnCount()>columnIndex){
			ColumnType type = screen.getColumnTypeManager().getColumnType(columnIndex);				
			PanelTableHeaderRenderer head = table.getHeaderRenderer(columnIndex);

			// Does color checking
			head.setTypeText(type.toString());
			
			// Identifier
			if(screen.getColumnTypeManager().getColumnType(columnIndex).equals(ColumnType.IDENTIFIER_LABEL)){
				if(screen.getColumnTypeManager().getCountOfType(ColumnType.IDENTIFIER_LABEL) > 1){
					// More than one identifier
					head.setTypeColor(Color.RED);
				} else {
					head.setTypeColor(Color.BLACK);
				}
			} 
			
			// Others than identifier
			else {
				if(screen.getColumnTypeManager().isChipNumberSetProperly(columnIndex)){
					// The ANNOTATION is always black, because the isChipNumberSetProperly return always 
					// true if column type is annotation
					head.setTypeColor(Color.BLACK);
					head.setTypeToolTipText("");
				} else {
					head.setTypeColor(Color.RED);
					if(screen.getColumnTypeManager().getChipCount() < screen.getColumnTypeManager().getCountOfType(type)){
						// User has selected more columns for this type than there is chips
						head.setTypeToolTipText(NOT_ENOUGH_CHIPS_TOOLTIP);
					} else {
						// The chip number is not set properly
						head.setTypeToolTipText(CHIP_NUMBER_NOT_SET_TOOLTIP);
					}
				}
			}
		}
	}

	public void actionPerformed(ActionEvent e) {

		// Listeners for first step
		if(screen.getCurrentStep() == ImportScreen.Step.FIRST){

			// Goto top button
			if (e.getSource() == toTopButton) {
				tableScroller.getVerticalScrollBar().setValue(
						tableScroller.getVerticalScrollBar().getMinimum());
			}

			// Goto bottom button
			else if (e.getSource() == toBottomButton) {
				tableScroller.getVerticalScrollBar().setValue(
						tableScroller.getVerticalScrollBar().getMaximum());
			}

			// Header button
			else if (e.getSource() == markHeaderButton) {
				if (markHeaderButton.isSelected()) {
					selectButton(markHeaderButton, markHeaderFooterButtons);				
				}
				updateHilighter();
			}

			// Footer button
			else if (e.getSource() == markFooterButton) {
				if (markFooterButton.isSelected()) {
					selectButton(markFooterButton, markHeaderFooterButtons);				
				}
				updateHilighter();
			}

			// Mark title button
			else if (e.getSource() == markTitleButton) {
				if (markTitleButton.isSelected()) {
					selectButton(markTitleButton, markHeaderFooterButtons);				
				}
				updateHilighter();
			} 
			
			// Reset button first step
			else if(e.getSource() == resetButtonFirstStep) {
				screen.getConversionModel().setColumnTitleLine(-1);
				screen.getConversionModel().setHeaderEndsRow(-1);
				screen.getConversionModel().setFooterBeginsRow(
						screen.getConversionModel().getLimitedRowCount());
			}
		}

		// Listeners for second step
		if(screen.getCurrentStep() == ImportScreen.Step.SECOND){

			// Deselect other buttons (step 2)
			if (markColumnsButtons != null && markColumnsButtons.containsKey(e.getSource())){
				selectButton((JToggleButton)e.getSource(), markColumnsButtons.keySet());
				updateHilighter();
			}
			
			// Reset button second step
			else if(e.getSource() == resetButtonSecondStep){
				screen.getColumnTypeManager().resetColumnTypes();
			}
		}
	}

	/**
	 * Checks if ImportScreen's mark header button is selected and 
	 * if it is then the table is in header marking mode
	 * 
	 * @return is table in header marking mode
	 */
	public boolean isInHeaderMarkingMode() {
		if (screen.getCurrentStep() == ImportScreen.Step.FIRST) {
			return this.markHeaderButton.isSelected();
		} else {
			return false;
		}
	}

	/**
	 * Checks if ImportScreen's mark footer button is selected and 
	 * if it is then the table is in footer marking mode
	 * 
	 * @return is table in footer marking mode
	 */
	public boolean isInFooterMarkingMode() {
		if (screen.getCurrentStep() == ImportScreen.Step.FIRST) {
			return this.markFooterButton.isSelected();
		} else {
			return false;
		}
	}

	/**
	 * Checks if ImportScreen's mark title button is selected and 
	 * if it is then the table is in title marking mode
	 * 
	 * @return is table in title marking mode
	 */
	public boolean isInTitleMarkingMode() {
		if (screen.getCurrentStep() == ImportScreen.Step.FIRST) {
			return this.markTitleButton.isSelected();
		} else {
			return false;
		}
	}
	
	/**
	 * Return a ColumnType depending on which button is pressed. 
	 * @return selected column type
	 */
	public ColumnType getSelectedColumnType(){
		if(markColumnsButtons == null){
			return null;
		}
		for(JToggleButton button : markColumnsButtons.keySet()){
			if(button.isSelected()){
				return markColumnsButtons.get(button);
			}
		}
		return null;
		
	}

	/**
	 * Sets spinners value to point to the row number of footer's start row.
	 * @param header
	 *            Footer start row (one bigger than table's value)
	 */
	public void setFooterStartRow(int footer) {
		// Set limits
		int header = Integer.parseInt(headerSpinner.getValue().toString());

		// Proper values are between header and total row count
		if (footer <= header) {
			footer = header + 1;
		} 
		else if (table != null && footer > table.getRowCount() + 1) {
			// The maximum value is row count plus one
			// The plus one means no footer at all
			footer = table.getRowCount() + 1;
		}

		this.footerSpinner.setValue(footer);
	}
	
	/**
	 * Resets the spinner values. Notice that calling this method will fire 
	 * an stateChanged event
	 *
	 */
	public void resetSpinners(){
		setHeaderEndRow(-1);
		setFooterStartRow(Integer.MAX_VALUE);
	}

	/**
	 * Sets spinners value to point to the row number of header's ending row.
	 * 
	 * @param header
	 *            Header end row (one bigger than table's value)
	 */
	public void setHeaderEndRow(int header) {
		int footer = Integer.parseInt(footerSpinner.getValue().toString());

		// Proper values are between 0 and footer
		if (header >= footer) {
			header = footer - 1;
		}
		else if (header < 0) {
			header = 0;
		}

		this.headerSpinner.setValue((header));
	}

	public void stateChanged(ChangeEvent e) {
		// Header spinner
		if (e.getSource() == headerSpinner) {
			int header = Integer.parseInt(headerSpinner.getValue().toString());
			
			// Do not change conversion model's values if it is correct already
			// Prevents duplicated events
			if(header -1 != screen.getConversionModel().getHeaderEnd()){
				screen.getConversionModel().setHeaderEndsRow(header -1);
			}
		}

		else if (e.getSource() == footerSpinner) {
			int footer = Integer.parseInt(footerSpinner.getValue().toString());
			
			// Do not change conversion model's values if it is correct already
			// Prevents duplicated events
			// Check also that the value is different than the default value Integer.MAX
			if(footer -1 != screen.getConversionModel().getFooterStart() && 
					footer != Integer.MAX_VALUE){
				screen.getConversionModel().setFooterBeginsRow(footer -1);
			}
		}
	}

	public void decimalSeparatorChanged(DecimalSeparatorChangedEvent e) {
		// Do nothing
	}

	public void delimiterChanged(DelimiterChangedEvent e) {
		// Do nothing
	}

	public void footerChanged(FooterChangedEvent e) {
		if(e.getNewValue() == -1){
			// User dragged mouse out of the table
			this.setFooterStartRow(table.getRowCount() + 1);
		} else {
			this.setFooterStartRow(e.getNewValue() +1);
		}
	}

	public void headerChanged(HeaderChangedEvent e) {
		this.setHeaderEndRow(e.getNewValue() +1);
	}

	public void titleRowChanged(TitleRowChangedEvent e) {
		// Do nothing
	}

	public void chipNumberChanged(ChipNumberChangedEvent event) {		
		updateChipCountComboBox(event.getColumnIndex());
		updateColumnTitleLabel(event.getColumnIndex());
	}

	public void columnTypeChanged(ColumnTypeChangedEvent event) {
		this.updateColumnTitleLabel(event.getColumnIndex());
	}

	public void chipCountChanged(ChipCountChangeEvent event) {
		updateAllChipCountComboBoxes();
		updateAllColumnTitleLabels();
	}

	public void columnTitlesChanged(ColumnTitlesChangedEvent e) {
		// Do nothing
	}

	public void inputFileChanged(InputFileChangedEvent e) {
		updateTitle();
	}
}
