package fi.csc.microarray.client.dataimport.tools;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXTaskPane;
import org.jdesktop.swingx.JXTaskPaneContainer;

import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.dataimport.ColumnTypePattern;
import fi.csc.microarray.client.dataimport.Delimiter;
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
 * Class for import tools (left side of the main split pane)
 * 
 * @author mkoski klemela
 *
 */
public class ToolsInternalFrame extends SimpleInternalFrame 
	implements ActionListener, CaretListener, ColumnTypeChangeListener, ConversionModelChangeListener{
	
	private static final Logger logger = Logger.getLogger(ToolsInternalFrame.class);

	//private static final String TYPE_RAW = "raw";
	//private static final String TYPE_NORM = "norm";
	
	private JXTaskPaneContainer firstStepOptionPanel;
	private JXTaskPaneContainer secondStepOptionPanel;
	
	private List<JRadioButton> delimRadioButtons;
	private JRadioButton customDelimRadioButton;
	private JTextField customDelimField;
	private JRadioButton dotAsDecimalSeparatorRadioButton;
	private JRadioButton commaAsDecimalSeparatorRadioButton;
	
	private ImportScreen screen;
	private ChipCountPanel chipCountPanel;

	private JButton fillTheRestButton;
	private JButton useCustomDelimButton;
	private JButton undoGuessButton;

	//private JTextField chipTypeField;
	//private String DEFAULT_CHIP_TYPE = "cDNA";
	//private ButtonGroup dataTypeGroup;
	//private JLabel chipTypeLabel;

	public ToolsInternalFrame(ImportScreen screen) {
		super("Tools");
		this.screen = screen;
	}
	
	/**
	 * Initializes the internal frame for the first step. Creates necessery 
	 * components and sets them to frame.
	 *
	 */
	public void initializeFirstStep(){		
		firstStepOptionPanel = new JXTaskPaneContainer();
		firstStepOptionPanel.add(this.createDelimSelectorPanel());
		firstStepOptionPanel.add(this.createDecimalSeparatorPanel());
		JScrollPane scroll = new JScrollPane(firstStepOptionPanel);
		scroll.setBorder(BorderFactory.createEmptyBorder());
		this.setContent(scroll);
	}
	
	/**
	 * Initializes the internal frame for the second step. Creates necesserly 
	 * components and sets them to frame.
	 *
	 */
	public void initializeSecondStep(){
		secondStepOptionPanel = new JXTaskPaneContainer();

//		secondStepOptionPanel.add(this.getDataTypePanel());
		secondStepOptionPanel.add(this.getChipCountPanel());		
		secondStepOptionPanel.add(this.createGuessTheRestPanel());		
		secondStepOptionPanel.add(this.createDataTrimmingPanel());
		JScrollPane scroll = new JScrollPane(secondStepOptionPanel);
		scroll.setBorder(BorderFactory.createEmptyBorder());
		this.setContent(scroll);
	}
	
	
	

	// F I R S T  S T E P ////////////////////////////////////////////////////////

	/**
	 * Creates panel for delimeter selection. 
	 * This panel is used on the first step
	 * 
	 * @return JPanel delimeter selection panel
	 */
	private JPanel createDelimSelectorPanel() {
		JXTaskPane delimPanel = new JXTaskPane();
		delimPanel.setLayout(new GridBagLayout());

		delimRadioButtons = new ArrayList<JRadioButton>();
		Map<String, JRadioButton>delimRadioButtonsByDelims = new HashMap<String, JRadioButton>();
		ButtonGroup delimGroup = new ButtonGroup();
		
		String[] selectorLabels = new String[Delimiter.values().length];
		String[] delimiters = new String[Delimiter.values().length];
		for(int i = 0; i < Delimiter.values().length; i++){
			selectorLabels[i] = Delimiter.values()[i].getName();
			delimiters[i] = Delimiter.values()[i].toString();
		}

		for (int i = 0; i < selectorLabels.length; i++) {
			JRadioButton selector = new JRadioButton(selectorLabels[i]);
			selector.setActionCommand(delimiters[i]);
			selector.addActionListener(this);
			selector.setOpaque(false);
			delimGroup.add(selector);
			delimRadioButtons.add(selector);
			delimRadioButtonsByDelims.put(delimiters[i], selector);
		}

		customDelimRadioButton = new JRadioButton("Other:");
		customDelimRadioButton.addActionListener(this);
		customDelimRadioButton.setOpaque(false);
		delimGroup.add(customDelimRadioButton);
		delimRadioButtons.add(customDelimRadioButton);

		delimRadioButtons.get(0).setSelected(true);

		customDelimField = new JTextField(3);
		customDelimField.addCaretListener(this);
		//customDelimField.setPreferredSize(new Dimension(25, 20));
		customDelimField.setMargin(new Insets(2, 2, 2, 2));
		
		useCustomDelimButton = new JButton("Use");
		useCustomDelimButton.addActionListener(this);
		useCustomDelimButton.setEnabled(false);

		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0; c.gridy = 0;
		c.gridwidth = 2;
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;

		for (int i = 0; i < delimRadioButtons.size()-1; i++) {
			delimPanel.add(delimRadioButtons.get(i), c);
			c.gridy++;			
		}

		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 0.0;
		c.fill = GridBagConstraints.NONE;
		delimPanel.add(customDelimRadioButton, c);
		c.gridx = 1;
		delimPanel.add(customDelimField, c);
		c.gridx++;
		delimPanel.add(useCustomDelimButton, c);

		delimPanel.setTitle("Column Delimiter");

		return delimPanel;
	}
	
	/**
	 * Creates panel for decimal separator selecting.
	 * This panel is used on the first step.
	 * 
	 * @return panel for decimal separator selecting
	 */
	private JPanel createDecimalSeparatorPanel() {
		JXTaskPane decimalSeparatorPanel = new JXTaskPane();
		decimalSeparatorPanel.setLayout(new GridBagLayout());

		ButtonGroup separatorGroup = new ButtonGroup();

		dotAsDecimalSeparatorRadioButton = new JRadioButton("Dot .");
		dotAsDecimalSeparatorRadioButton.setActionCommand(".");
		dotAsDecimalSeparatorRadioButton.addActionListener(this);
		dotAsDecimalSeparatorRadioButton.setOpaque(false);
		separatorGroup.add(dotAsDecimalSeparatorRadioButton);

		commaAsDecimalSeparatorRadioButton = new JRadioButton("Comma ,");
		commaAsDecimalSeparatorRadioButton.setActionCommand(",");
		commaAsDecimalSeparatorRadioButton.addActionListener(this);
		commaAsDecimalSeparatorRadioButton.setOpaque(false);
		separatorGroup.add(commaAsDecimalSeparatorRadioButton);

		dotAsDecimalSeparatorRadioButton.setSelected(true);
		screen.getConversionModel().setDecimalSeparator('.');

		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0; c.gridy = 0;
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;

		decimalSeparatorPanel.add(dotAsDecimalSeparatorRadioButton, c);
		c.gridy++;
		decimalSeparatorPanel.add(commaAsDecimalSeparatorRadioButton, c);

		decimalSeparatorPanel.setTitle("Decimal Separator");

		return decimalSeparatorPanel;
	}
	
	
	/**
	 * Sets delimeter for values of a row. Checks for collision with decimalSeparator.
	 * @param delim
	 */
	private void setDelim(String delim) {
		screen.getConversionModel().setDelim(Delimiter.stringToDelim(delim));
		if (screen.getConversionModel().getDecimalSeparator() == 
			screen.getConversionModel().getDelim().toString().charAt(0) &&
			screen.getConversionModel().getDelim().toString().length() == 1) {
			
			JOptionPane.showMessageDialog(null,
					new JLabel("<html>You can't use the same character<br>as both the column delimiter and the<br>decimal separator!</html>"),
					"Error",
					JOptionPane.ERROR_MESSAGE);
			customDelimField.setForeground(Color.RED);
			useCustomDelimButton.setEnabled(false);
			
		} else {
		
			this.updateDelimeterPanel();
				
			try {
				// Test pattern validity:
				Pattern.compile(delim);
				customDelimField.setForeground(Color.BLACK);
				if(Delimiter.isCustom(delim)){
					useCustomDelimButton.setEnabled(true);
				}
			} catch (PatternSyntaxException pse) {
				customDelimField.setForeground(Color.RED);
				if(Delimiter.isCustom(delim)){
					useCustomDelimButton.setEnabled(false);
				}
			}
		}
	}

	private void setDecimalSeparator(char separator) {
		screen.getConversionModel().setDecimalSeparator(separator);
		if (screen.getConversionModel().getDecimalSeparator() == 
			screen.getConversionModel().getDelim().toString().charAt(0) &&
			screen.getConversionModel().getDelim().toString().length() == 1) {
			
			JOptionPane.showMessageDialog(null,
					new JLabel("<html>You can't use the same character<br>as both the column delimiter and the<br>decimal separator!</html>"),
					"Error",
					JOptionPane.ERROR_MESSAGE);
			customDelimField.setForeground(Color.RED);
			useCustomDelimButton.setEnabled(false);
			
		} else {
		
			this.updateDelimeterPanel();
			customDelimField.setForeground(Color.BLACK);
			useCustomDelimButton.setEnabled(true);				
		}
	}

	// S E C O N D  S T E P ///////////////////////////////////////////////////////
	
	/**
	 * Panel that should change the column headers to enable import of normalised
	 * datas in the future. For
	 * version 1.4.0 this functionality isn't ready yet and this panel is just hidden. 
	 * 
	 * @return
	 */
//	private JXTaskPane getDataTypePanel() {
//
//		JXTaskPane dataTypePanel = new JXTaskPane();
//		dataTypePanel.setTitle("Data type");
//		dataTypePanel.setLayout(new BorderLayout());
//
//
//		dataTypePanel.setLayout(new GridBagLayout());
//
//		dataTypeGroup = new ButtonGroup();		
//		JRadioButton rawButton = new JRadioButton("Raw data");
//		JRadioButton normButton = new JRadioButton("Normalised data");
//		rawButton.setOpaque(false);
//		normButton.setOpaque(false);
//		rawButton.setActionCommand(TYPE_RAW);
//		normButton.setActionCommand(TYPE_NORM);
//		dataTypeGroup.add(rawButton);
//		dataTypeGroup.add(normButton);
//
//		rawButton.addActionListener(new ActionListener(){
//			public void actionPerformed(ActionEvent e) {
//				setNormalised(false);
//			}	
//		});
//
//		normButton.addActionListener(new ActionListener(){
//			public void actionPerformed(ActionEvent e) {
//				setNormalised(true);
//			}	
//		});
//
//		chipTypeLabel = new JLabel("Chiptype:");
//		chipTypeField = new JTextField(DEFAULT_CHIP_TYPE );
//		chipTypeField.setEnabled(false);
//		chipTypeLabel.setEnabled(false);
//		rawButton.setSelected(true);	
//
//		//chipTypeField.setMargin(new Insets(2, 2, 2, 2));
//
//		GridBagConstraints c = new GridBagConstraints();
//		c.gridx = 0; c.gridy = 0;
//		c.gridwidth = 1;
//		c.anchor = GridBagConstraints.WEST;
//		c.weightx = 1.0;
//		c.fill = GridBagConstraints.HORIZONTAL;
//
//		dataTypePanel.add(rawButton, c);
//		c.gridy++;
//		dataTypePanel.add(normButton, c);
//		c.gridy++;
//		c.insets.top += 5;
//		dataTypePanel.add(chipTypeLabel, c);
//		c.insets.top -= 5;
//		c.gridy++;
//		dataTypePanel.add(chipTypeField, c);
//
//		return dataTypePanel;
//
//	}
//
//	private void setNormalised(boolean isNormalised){
//		chipTypeField.setEnabled(isNormalised);
//		chipTypeLabel.setEnabled(isNormalised);
//		
//		screen.getConversionModel().setNormalised(isNormalised);
//	}
	
	public ChipCountPanel getChipCountPanel() {
		if(chipCountPanel == null){
			chipCountPanel = new ChipCountPanel(screen);
			chipCountPanel.setTitle("Chip counts");
			chipCountPanel.setToolTipText("Shows the count of selected columns of each type");
		}
		return chipCountPanel;
	}
	
	/**
	 * Find-and-replace functionality.
	 */
	private DataTrimmingPanel createDataTrimmingPanel() {
		DataTrimmingPanel dataTrimmingPanel = new DataTrimmingPanel(screen);
		screen.getConversionModel().addConversionChangeListener(dataTrimmingPanel);
		return dataTrimmingPanel;
	}

	
	private Component createGuessTheRestPanel() {
		JXTaskPane guessTheRestPanel = new JXTaskPane();
		guessTheRestPanel.setLayout(new BorderLayout());
		fillTheRestButton = new JButton("Complete the rest");
		guessTheRestPanel.add(fillTheRestButton, BorderLayout.WEST);
		fillTheRestButton.addActionListener(this);
		guessTheRestPanel.setTitle("Complete with pattern");
		undoGuessButton = new JButton("Undo");
		undoGuessButton.addActionListener(this);
		undoGuessButton.setEnabled(false);
		guessTheRestPanel.add(undoGuessButton, BorderLayout.EAST);
		return guessTheRestPanel;
	}
	
	
	// L I S T E N E R S //////////////////////////////////////////////////////
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source instanceof JRadioButton) {
			JRadioButton button = (JRadioButton) source;
			
			// Delimeter selection radio button
			if (delimRadioButtons.contains(button)) {
				
				// Predefined delimeter
				if (button != customDelimRadioButton) {					
					this.setDelim(button.getActionCommand());			
					screen.updateTable(false);
				} 
				
				// Custom delimeter
				else {
					String customDelim = customDelimField.getText();
					if (customDelim != null && customDelim.length() > 0) {
						this.setDelim(customDelim);
					}
				}
			} 
			
			// Decimal separator selection
			else if (button == dotAsDecimalSeparatorRadioButton ||
					button == commaAsDecimalSeparatorRadioButton) {
				screen.getConversionModel().setDecimalSeparator(button.getActionCommand().charAt(0));
			}
		}
		
		if(source == fillTheRestButton){
			ColumnTypePattern pattern = ColumnTypePattern.createColumnTypePatternFromAllColumns(screen.getColumnTypeManager().getColumns());
			screen.getColumnTypeManager().selectColumnsFromPattern(pattern, true, screen.getConversionModel().getCleanColumnTitles());
			undoGuessButton.setEnabled(true);
			screen.getTableFrame().getTable().repaint();
		}
		
		if(source == undoGuessButton){
			screen.getColumnTypeManager().undoPatternSelection();
			undoGuessButton.setEnabled(false);
		}
		
		if(source == useCustomDelimButton){
			screen.updateTable(false);
			
			// Disable button after data parsing and enable it again when the delimiter 
			// changes
			useCustomDelimButton.setEnabled(false);
		}
	}

	public void caretUpdate(CaretEvent e) {
		Object source = e.getSource();
		
		// Custom delimeter field
		if (source == customDelimField) {
			if (customDelimRadioButton.isSelected()) {
				String customDelim = customDelimField.getText();
				if (customDelim != null && customDelim.length() > 0) {
					this.setDelim(customDelim);
				} else {
					customDelimField.setForeground(Color.BLACK);
				}	
			}
		}
	}

	public void updateDelimeterPanel() {
		char decimalSeparator = screen.getConversionModel().getDecimalSeparator();
		Delimiter delim = screen.getConversionModel().getDelim();
		if(decimalSeparator == '.'){
			dotAsDecimalSeparatorRadioButton.setSelected(true);
			commaAsDecimalSeparatorRadioButton.setSelected(false);
		} else {
			dotAsDecimalSeparatorRadioButton.setSelected(false);
			commaAsDecimalSeparatorRadioButton.setSelected(true);
		}
		
		commaAsDecimalSeparatorRadioButton.setEnabled(delim != Delimiter.COMMA);
		
		//To show right radiobutton selected after fileAnalyzer
		for(JRadioButton button : this.delimRadioButtons ){
			button.setSelected(button.getActionCommand().equals(delim.toString()));
		}
	}
	
	public void columnTypeChanged(ColumnTypeChangedEvent event) {
		
		logger.debug("columnTypeChanged");
		this.getChipCountPanel().updateAllKeyColumnCounters();		
	}

	public void chipNumberChanged(ChipNumberChangedEvent event) {
		this.getChipCountPanel().updateAllKeyColumnCounters();
	}

	public void chipCountChanged(ChipCountChangeEvent event) {
		// Do nothing
	}

	/**
	 * Updates GUI if decimal separator changed
	 */
	public void decimalSeparatorChanged(DecimalSeparatorChangedEvent e) {
		this.setDecimalSeparator(e.getNewValue().toString().charAt(0));
	}

	/**
	 * Updates GUI if delimiter changed
	 */
	public void delimiterChanged(DelimiterChangedEvent e) {
		this.setDelim(e.getNewValue().toString());
	}

	public void footerChanged(FooterChangedEvent e) {
		// Do nothing
	}

	public void headerChanged(HeaderChangedEvent e) {
		// Do nothing
	}

	public void titleRowChanged(TitleRowChangedEvent e) {
		// Do nothing
	}

	public void columnTitlesChanged(ColumnTitlesChangedEvent e) {
		// Do nothing
	}

	public void inputFileChanged(InputFileChangedEvent e) {
		// Do nothing
	}
}
