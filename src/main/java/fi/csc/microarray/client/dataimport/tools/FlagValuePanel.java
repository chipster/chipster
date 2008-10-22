package fi.csc.microarray.client.dataimport.tools;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JComboBox;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.text.InternationalFormatter;

import org.jdesktop.swingx.JXTaskPane;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.dataimport.ColumnType;
import fi.csc.microarray.client.dataimport.ColumnTypeManager;
import fi.csc.microarray.client.dataimport.DataColumn;
import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.trimmer.ConditionalStringReplace;
import fi.csc.microarray.client.dataimport.trimmer.DataTrimmingOperation;
import fi.csc.microarray.client.dataimport.trimmer.NormalStringReplace;
import fi.csc.microarray.client.dataimport.trimmer.NumberToStringReplace;

public class FlagValuePanel extends JXTaskPane implements ActionListener,
        CaretListener, ChangeListener {

	private JLabel scaleLabel;
	private JComboBox scaleSelector;
    private JLabel[] flagNameLabels;
    
    private JPanel topFlagBoundariesPanel;
    private JPanel middleFlagBoundariesPanel;
    private JPanel bottomFlagBoundariesPanel;
    private JTextField topFlagValuesField;
    private JTextField middleFlagValuesField;
    private JTextField bottomFlagValuesField;
    private JComboBox topFlagOrderSelector;
    private JLabel bottomFlagOrderLabel;
    private JSpinner topFlagBoundarySpinner;
    private JFormattedTextField topFlagBoundaryField; 
    private JLabel middleFlagLowerBoundaryLabel;
    private JLabel middleBoundaryLabel;
    private JLabel middleFlagUpperBoundaryLabel;
    private JSpinner bottomFlagBoundarySpinner;
    private JFormattedTextField bottomFlagBoundaryField; 
    
    private String scale;
    private boolean showingDiscreteScaleExample;
    private double lowerBoundary;
    private double upperBoundary;
    private JSpinner lowerBoundarySpinner;
    private JSpinner upperBoundarySpinner;
    private boolean boundariesAreValid;
    private Map<JFormattedTextField, Integer> extraTokensForBoundaryFields;
    private Map<JFormattedTextField, JSpinner> spinnersForBoundaryFields;
    
    private DecimalFormat decimalFormat;
    
    private String topFlagDefault;
    private String middleFlagDefault;
    private String bottomFlagDefault;
    private double topFlagDefaultBoundary;
    private double bottomFlagDefaultBoundary;
    
    private ImportScreen screen;
    
    public static final Color TOP_FLAG_COLOR = new Color(210, 255, 200); // green
    public static final Color MIDDLE_FLAG_COLOR = new Color(255, 255, 180); // yellow
    public static final Color BOTTOM_FLAG_COLOR = new Color(255, 210, 200); // red
    
    public static final String DISABLED_SCALE;
    public static final String CONTINUOUS_SCALE;
    public static final String DISCRETE_SCALE;
    
    private static final String LOWER_BOUNDARY_LABEL;
    private static final String UPPER_BOUNDARY_LABEL;
    private static final String MIDDLE_BOUNDARY_LABEL;
    
    private static final int ORDER_SELECTOR_WIDTH;
    
    static {
    	LOWER_BOUNDARY_LABEL = "x <";
    	UPPER_BOUNDARY_LABEL = "x \u2265";
    	MIDDLE_BOUNDARY_LABEL = "\u2264 x <";
    	
    	ORDER_SELECTOR_WIDTH = 50;
    	
    	DISABLED_SCALE = "Disabled";
    	CONTINUOUS_SCALE = "Continuous";
    	DISCRETE_SCALE = "Discrete";
    }
    
    public FlagValuePanel(ImportScreen screen) {
        this("Present P:", "Marginal M", "Absent A",
                "P",       "M",       "A",
                     0.0,       0.0, screen);
        //was 0.5, 0.48
    }
    
    public FlagValuePanel(String topFlagName, String middleFlagName, String bottomFlagName,
            String topFlagDefault, String middleFlagDefault, String bottomFlagDefault,
            double topFlagDefaultBoundary, double bottomFlagDefaultBoundary,
            ImportScreen screen) {
        
        this.setLayout(new GridBagLayout());
        
        // To remember the defaults:
        this.topFlagDefault = topFlagDefault;
        this.middleFlagDefault = middleFlagDefault;
        this.bottomFlagDefault = bottomFlagDefault;
        this.topFlagDefaultBoundary = topFlagDefaultBoundary;
        this.bottomFlagDefaultBoundary = bottomFlagDefaultBoundary;

        this.screen = screen;
        
        Dimension selectorSize = new Dimension(ORDER_SELECTOR_WIDTH, 20);
        Dimension boundarySpinnerSize = new Dimension(63, 20);
        Dimension valueFieldSize = new Dimension(113, 22);
        
        scaleLabel = new JLabel("Scale:");
        scaleSelector = new JComboBox(
        		new String[] { DISABLED_SCALE, CONTINUOUS_SCALE, DISCRETE_SCALE } );
        scaleSelector.addActionListener(this);
        
        flagNameLabels = new JLabel[3];
        flagNameLabels[0] = new JLabel(topFlagName);
        flagNameLabels[1] = new JLabel(middleFlagName);
        flagNameLabels[2] = new JLabel(bottomFlagName);            
        
        topFlagOrderSelector = new JComboBox(
                new String[] { LOWER_BOUNDARY_LABEL, UPPER_BOUNDARY_LABEL } );
        topFlagOrderSelector.addActionListener(this);
        topFlagOrderSelector.setPreferredSize(selectorSize);
        bottomFlagOrderLabel = new JLabel();
        
        extraTokensForBoundaryFields = new HashMap<JFormattedTextField, Integer>();
        spinnersForBoundaryFields = new HashMap<JFormattedTextField, JSpinner>();
        
        topFlagBoundarySpinner = new JSpinner(new SpinnerNumberModel(
                0, -Double.MAX_VALUE, Double.MAX_VALUE, 0.01));
        topFlagBoundarySpinner.addChangeListener(this);
        topFlagBoundarySpinner.setPreferredSize(boundarySpinnerSize);
        topFlagBoundaryField = ((JSpinner.DefaultEditor) topFlagBoundarySpinner.getEditor()).getTextField();
        topFlagBoundaryField.addCaretListener(this);
        ((InternationalFormatter) topFlagBoundaryField.getFormatter()).setFormat(this.getDecimalFormat());
        spinnersForBoundaryFields.put(topFlagBoundaryField, topFlagBoundarySpinner);
        
        bottomFlagBoundarySpinner = new JSpinner(new SpinnerNumberModel(
                0, -Double.MAX_VALUE, Double.MAX_VALUE, 0.01));
        bottomFlagBoundarySpinner.addChangeListener(this);
        bottomFlagBoundarySpinner.setPreferredSize(boundarySpinnerSize);
        bottomFlagBoundaryField = ((JSpinner.DefaultEditor) bottomFlagBoundarySpinner.getEditor()).getTextField();
        bottomFlagBoundaryField.addCaretListener(this);
        ((InternationalFormatter) bottomFlagBoundaryField.getFormatter()).setFormat(this.getDecimalFormat());
        spinnersForBoundaryFields.put(bottomFlagBoundaryField, bottomFlagBoundarySpinner);
        
        topFlagBoundariesPanel = new JPanel(new GridBagLayout());
        middleFlagBoundariesPanel = new JPanel(new GridBagLayout());
        bottomFlagBoundariesPanel = new JPanel(new GridBagLayout());
        
        topFlagBoundariesPanel.setMinimumSize(valueFieldSize);
        middleFlagBoundariesPanel.setMinimumSize(valueFieldSize);
        bottomFlagBoundariesPanel.setMinimumSize(valueFieldSize);
        middleFlagBoundariesPanel.setOpaque(false);
        bottomFlagBoundariesPanel.setOpaque(false);  
        
        middleFlagLowerBoundaryLabel = new JLabel();
        middleFlagUpperBoundaryLabel = new JLabel();
        middleBoundaryLabel = new JLabel(MIDDLE_BOUNDARY_LABEL);
        middleBoundaryLabel.setOpaque(false);
        
        topFlagValuesField = new JTextField();
        middleFlagValuesField = new JTextField();
        bottomFlagValuesField = new JTextField();
        
        // At first, these will show a short help message:
        this.setShowingDiscreteScaleExample(true);

        // start listening only after components are properly initialised 
        bottomFlagValuesField.addCaretListener(this);
        middleFlagValuesField.addCaretListener(this);
        topFlagValuesField.addCaretListener(this);

        this.showScale(DISABLED_SCALE);
        this.applyDefaultBoundaries();        
        this.setFlagModificationEnabled(false);
        
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1;        
        c.weighty = 0;
        c.gridx = 0; 
        c.gridy = 0;
        c.anchor = GridBagConstraints.NORTHWEST;

        this.add(scaleLabel, c);
        c.gridx++;    
        this.add(scaleSelector, c);
        
        c.gridx = 0; 
        c.gridy = 1;
        c.gridwidth = 2;
        this.add(flagNameLabels[0], c);
        c.gridy += 2;
        this.add(flagNameLabels[1], c);
        c.gridy += 2;
        this.add(flagNameLabels[2], c);
               
        c.gridx = 0; 
        c.gridy = 2;
        c.weighty = 1;
        c.insets.bottom = 10;
        this.add(topFlagBoundariesPanel, c);
        c.gridy += 2;
        this.add(middleFlagBoundariesPanel, c);
        c.gridy += 2;
        this.add(bottomFlagBoundariesPanel, c);
        
        this.updateMiddleFlagBoundaries();
    }
    
    public void setScale(String scaleId) {
    	scaleSelector.setSelectedItem(scaleId);
    	updateFlagTrimmer();
    	// will also result in calling showScale
    }
    
    private void showScale(String scaleId) {
    	if (scaleId != null && (this.scale == null || !this.scale.equals(scaleId))) {
    		    		    		
    		this.scale = scaleId;
	    	topFlagBoundariesPanel.removeAll();
	    	middleFlagBoundariesPanel.removeAll();
	    	bottomFlagBoundariesPanel.removeAll();
	    	
	        GridBagConstraints c = new GridBagConstraints();
	        
	        c.fill = GridBagConstraints.HORIZONTAL;
	        
	    	if (scaleId.equals(CONTINUOUS_SCALE)) {
	    		this.setFlagModificationEnabled(true);
	    		c.gridx = 0;
	    		c.gridy = 0;
	    		c.weightx = 1.0;
		        topFlagBoundariesPanel.add(topFlagOrderSelector, c);
		        		        
		        bottomFlagBoundariesPanel.add(bottomFlagOrderLabel, c);
		        
		        c.gridx++;
		        c.weightx = 0;
		        topFlagBoundariesPanel.add(topFlagBoundarySpinner, c);
		        		    
		        bottomFlagBoundariesPanel.add(bottomFlagBoundarySpinner, c);
		        
		        c.gridx = 0; 
		        c.gridy = 0;
		        
		        middleFlagBoundariesPanel.add(middleFlagLowerBoundaryLabel, c);
		        c.gridx = 1;
		        middleFlagBoundariesPanel.add(middleBoundaryLabel, c);
		        c.gridx = 2;
		        middleFlagBoundariesPanel.add(middleFlagUpperBoundaryLabel, c);
	    	} else if (scaleId.equals(DISCRETE_SCALE)){
	    		this.setFlagModificationEnabled(true);
	    		c.weightx = 1.0;
		        topFlagBoundariesPanel.add(topFlagValuesField, c);
		        middleFlagBoundariesPanel.add(middleFlagValuesField, c);
		        bottomFlagBoundariesPanel.add(bottomFlagValuesField, c);
	    	} else if(scaleId.equals(DISABLED_SCALE)){
    			this.setFlagModificationEnabled(false);
    		}
    	}
    	this.revalidate();
    	this.repaint();
    }
    
    private void setShowingDiscreteScaleExample(boolean showing) {
    	if (showing) {
            showingDiscreteScaleExample = true;
            topFlagValuesField.setText("use ; as separator");
            middleFlagValuesField.setText("0;2;-1");
            bottomFlagValuesField.setText("");
    	} else {
			showingDiscreteScaleExample = false;
			topFlagValuesField.setText("");
			middleFlagValuesField.setText("");
			bottomFlagValuesField.setText("");
			topFlagValuesField.setForeground(Color.black);
			middleFlagValuesField.setForeground(Color.black);
			bottomFlagValuesField.setForeground(Color.black);
    	}
    }
    
    private DecimalFormat getDecimalFormat() {
        if (decimalFormat == null) {
            decimalFormat = new DecimalFormat("0.00##");
            DecimalFormatSymbols symbols = decimalFormat.getDecimalFormatSymbols();
            symbols.setDecimalSeparator('.');
            decimalFormat.setDecimalFormatSymbols(symbols);
        }
        return decimalFormat;
    }
    
    private void applyDefaultBoundaries() {
        topFlagBoundarySpinner.setValue(topFlagDefaultBoundary);
        bottomFlagBoundarySpinner.setValue(bottomFlagDefaultBoundary);
        boundariesAreValid = true;
        
        if (topFlagDefaultBoundary >= bottomFlagDefaultBoundary) {
        	this.setUpperBoundary(topFlagDefaultBoundary);
        	this.setLowerBoundary(bottomFlagDefaultBoundary);
            upperBoundarySpinner = topFlagBoundarySpinner;
            lowerBoundarySpinner = bottomFlagBoundarySpinner;
            topFlagOrderSelector.setSelectedItem(UPPER_BOUNDARY_LABEL);
            // fires also the ActionEvent which will result in calling setOrder
        } else {
        	this.setUpperBoundary(bottomFlagDefaultBoundary);
        	this.setLowerBoundary(topFlagDefaultBoundary);
            upperBoundarySpinner = bottomFlagBoundarySpinner;
            lowerBoundarySpinner = topFlagBoundarySpinner;
            topFlagOrderSelector.setSelectedItem(LOWER_BOUNDARY_LABEL);
        }
    }
    
    private void setOrder(String topFlagBoundaryLabel) {
        double upperBoundary = (Double) upperBoundarySpinner.getValue();
        double lowerBoundary = (Double) lowerBoundarySpinner.getValue();
        
        if (topFlagBoundaryLabel.equals(UPPER_BOUNDARY_LABEL)) {
            upperBoundarySpinner = topFlagBoundarySpinner;
            lowerBoundarySpinner = bottomFlagBoundarySpinner;
            bottomFlagOrderLabel.setText(LOWER_BOUNDARY_LABEL);
        } else if (topFlagBoundaryLabel.equals(LOWER_BOUNDARY_LABEL)) {
            upperBoundarySpinner = bottomFlagBoundarySpinner;
            lowerBoundarySpinner = topFlagBoundarySpinner;
            bottomFlagOrderLabel.setText(UPPER_BOUNDARY_LABEL);
        }
        if (boundariesAreValid) {
            upperBoundarySpinner.setValue(upperBoundary);
            lowerBoundarySpinner.setValue(lowerBoundary);
        } else {
            this.checkAndShowBoundaryValidity();
        }
        
        this.updateMiddleFlagBoundaries();
    }
    
    private void updateMiddleFlagBoundaries() {
        if (lowerBoundarySpinner != null && upperBoundarySpinner != null) {
            String lowerBoundaryText = this.getDecimalFormat().format(lowerBoundary);
            String upperBoundaryText = this.getDecimalFormat().format(upperBoundary);
            if (boundariesAreValid && !lowerBoundaryText.equals(upperBoundaryText)) {
                middleFlagLowerBoundaryLabel.setText(lowerBoundaryText);
                middleBoundaryLabel.setText(MIDDLE_BOUNDARY_LABEL);
                middleFlagUpperBoundaryLabel.setText(upperBoundaryText);
            } else {
                middleFlagLowerBoundaryLabel.setText(" ");
                middleBoundaryLabel.setText(" ");
                middleFlagUpperBoundaryLabel.setText(" ");
            }	
        }
    }
      
    private void checkAndShowBoundaryValidity() {
    	if (upperBoundary < lowerBoundary) {
    		if (boundariesAreValid) {
    			((JSpinner.DefaultEditor) upperBoundarySpinner.getEditor()).getTextField().setForeground(Color.RED);
    			((JSpinner.DefaultEditor) lowerBoundarySpinner.getEditor()).getTextField().setForeground(Color.RED);
    			boundariesAreValid = false;
    		}
    	} else {
    		if (!boundariesAreValid) {
    			((JSpinner.DefaultEditor) upperBoundarySpinner.getEditor()).getTextField().setForeground(Color.BLACK);
    			((JSpinner.DefaultEditor) lowerBoundarySpinner.getEditor()).getTextField().setForeground(Color.BLACK);
    			boundariesAreValid = true;
    		}
    	}
    }
    
    //TODO This method propably useless
    private void handleContentChangeFor(JFormattedTextField field) {
		String text = field.getText();
		try {
			// Might throw NumberFormatException:
			double value = Double.parseDouble(text);
			
			// From now on, we know that the input is a valid double:
			field.setForeground(Color.black);
			
			int textLength = text.length();
			int extraTokens = 0;
			if (textLength > 6) {
				extraTokens = (textLength-6);
				if (extraTokens > 10) {
					// extraTokens = 10;
				}
			}
			extraTokensForBoundaryFields.put(field, extraTokens);
			for (int archivedTokenCount : extraTokensForBoundaryFields.values()) {
				if (archivedTokenCount > extraTokens) {
					extraTokens = archivedTokenCount;
				}
			}
			
			if (spinnersForBoundaryFields.get(field) == lowerBoundarySpinner) {
				this.setLowerBoundary(value);
			} else {
				this.setUpperBoundary(value);
			}
			this.updateMiddleFlagBoundaries();
			this.checkAndShowBoundaryValidity();
			
			this.revalidate();
		} catch (NumberFormatException nfe) {
			field.setForeground(Color.red);
		}
    }
    
    private void setUpperBoundary(double value) {
    	upperBoundary = value;
    	updateFlagTrimmer();
    	screen.getTableFrame().getTable().repaint();
    }
    
    private void setLowerBoundary(double value) {
    	lowerBoundary = value;
    	updateFlagTrimmer();
    	screen.getTableFrame().getTable().repaint();
    }
    
    public String getTopFlag() {
        return topFlagDefault;
    }
    
    public String getMiddleFlag() {
        return middleFlagDefault;
    }
    
    public String getBottomFlag() {
        return bottomFlagDefault;
    }
    
    public double getTopFlagBoundary() {
        return (Double) topFlagBoundarySpinner.getValue();
    }
    
    public double getBottomFlagBoundary() {
        return (Double) bottomFlagBoundarySpinner.getValue();
    }
    
    public boolean topFlagHasBiggestValues() {
        return ((String) topFlagOrderSelector.getSelectedItem()).equals(UPPER_BOUNDARY_LABEL);
    }
    
    public void setTopFlagBoundary(double value) {
    	if (this.topFlagHasBiggestValues()) {
    		this.setUpperBoundary(value);
    		upperBoundarySpinner.setValue(value);
    	} else {
    		this.setLowerBoundary(value);
    		lowerBoundarySpinner.setValue(value);
    	}
    }
    
    public void setBottomFlagBoundary(double value) {
    	if (this.topFlagHasBiggestValues()) {
    		this.setLowerBoundary(value);
    		lowerBoundarySpinner.setValue(value);
    	} else {
    		this.setUpperBoundary(value);
    		upperBoundarySpinner.setValue(value);
    	}
    }
    
    public boolean isFlagModificationEnabled(){
    	//All components should have the same state, so any of them can be used
    	return flagNameLabels[0].isEnabled();
    }
    
    public void setEnabled(boolean enabled) {
    	scaleLabel.setEnabled(enabled);
    	scaleSelector.setEnabled(enabled);
    	if(enabled){
    		setScale(DISCRETE_SCALE);
    	} else {
    		setScale(DISABLED_SCALE);
    	}
    }
    
    public void setFlagModificationEnabled(boolean enabled){
    	for (int i = 0; i < flagNameLabels.length; i++) {
    		flagNameLabels[i].setEnabled(enabled);
    	}
    	topFlagValuesField.setEnabled(enabled);
    	middleFlagValuesField.setEnabled(enabled);
    	bottomFlagValuesField.setEnabled(enabled);
    	this.setBoundaryAdjustmentEnabled(enabled);
}
    
    public void setBoundaryAdjustmentEnabled(boolean enabled) {
        topFlagOrderSelector.setEnabled(enabled);
        topFlagBoundarySpinner.setEnabled(enabled);
        middleFlagLowerBoundaryLabel.setEnabled(enabled);
        middleBoundaryLabel.setEnabled(enabled);
        middleFlagUpperBoundaryLabel.setEnabled(enabled);
        bottomFlagOrderLabel.setEnabled(enabled);
        bottomFlagBoundarySpinner.setEnabled(enabled);	
    }
    
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        if (source == topFlagOrderSelector) {
            this.setOrder((String) topFlagOrderSelector.getSelectedItem());
        } else if (source == scaleSelector) {
        	this.showScale((String) scaleSelector.getSelectedItem());
        }
    }
    
    public void caretUpdate(CaretEvent e) {
    	Object source = e.getSource();
    	if (source instanceof JFormattedTextField) {
    		JFormattedTextField field = (JFormattedTextField) source;
    		this.handleContentChangeFor(field);
    	} else if (source instanceof JTextField) {
    		if (source == topFlagValuesField ||
    				source == middleFlagValuesField ||
    				source == bottomFlagValuesField) {
    			if (showingDiscreteScaleExample) {
    				this.setShowingDiscreteScaleExample(false);
    			}
    		}
    		updateFlagTrimmer();
    	}
    }

    public void stateChanged(ChangeEvent e) {
    	Object source = e.getSource();
    	if (source instanceof JSpinner) {
    		if (source == upperBoundarySpinner) {
    			this.setUpperBoundary((Double) upperBoundarySpinner.getValue());
    		} else if (source == lowerBoundarySpinner) {
    			this.setLowerBoundary((Double) lowerBoundarySpinner.getValue());
    		}
    		
    		if (upperBoundarySpinner != null && lowerBoundarySpinner != null) {
    			this.checkAndShowBoundaryValidity();
    		}
    		this.updateMiddleFlagBoundaries();
    		screen.getTableFrame().getTable().repaint();
    	}
    }
    
    public void updateFlagTrimmer(){
    	ColumnTypeManager manager = screen.getColumnTypeManager();
    	List<DataTrimmingOperation> operations = new ArrayList<DataTrimmingOperation>();
    	for(DataColumn column : manager.getColumns()){
    		if(column.getColumnType() == ColumnType.FLAG_LABEL){
    			if(scale.equals(CONTINUOUS_SCALE)){
    				int columnIndex = column.getColumnIndex();
    				if(topFlagHasBiggestValues()){
    					operations.add(new ConditionalStringReplace(Double.NEGATIVE_INFINITY, lowerBoundary, false, false, getBottomFlag(), columnIndex));
    					operations.add(new ConditionalStringReplace(lowerBoundary, upperBoundary, true, false, getMiddleFlag(), columnIndex));
    					operations.add(new ConditionalStringReplace(upperBoundary, Double.POSITIVE_INFINITY, true, false, getTopFlag(), columnIndex));
    				} else {
    					operations.add(new ConditionalStringReplace(Double.NEGATIVE_INFINITY, lowerBoundary, false, false, getTopFlag(), columnIndex));
    					operations.add(new ConditionalStringReplace(lowerBoundary, upperBoundary, true, false, getMiddleFlag(), columnIndex));
    					operations.add(new ConditionalStringReplace(upperBoundary, Double.POSITIVE_INFINITY, true, false, getBottomFlag(), columnIndex));
    				}
    			} else if(this.scale.equals(DISCRETE_SCALE)){
    				int columnIndex = column.getColumnIndex();
    				
    				if(topFlagValuesField.getText() != null && topFlagValuesField.getText().length() > 0){
    					operations.addAll(FlagValuePanel.parseDiscreteString(topFlagValuesField.getText(), getTopFlag(), columnIndex));
    				}
    				if(middleFlagValuesField.getText() != null && middleFlagValuesField.getText().length() > 0){
    					operations.addAll(FlagValuePanel.parseDiscreteString(middleFlagValuesField.getText(), getMiddleFlag(), columnIndex));
    				}
    				if(bottomFlagValuesField.getText() != null && bottomFlagValuesField.getText().length() > 0){
    					operations.addAll(FlagValuePanel.parseDiscreteString(bottomFlagValuesField.getText(), getBottomFlag(), columnIndex));
    				}
    			} 
    		}
    	}
    	screen.getFlagTrimmer().clear();
    	screen.getFlagTrimmer().pushOperations(operations);
    	screen.getTableFrame().getTable().repaint();
    }

	public static void main(String[] args) {
        // set the UI defaults
        UIDefaults defaults = UIManager.getDefaults();
        defaults.putAll(VisualConstants.getUIDefaults());
        
        FlagValuePanel panel = new FlagValuePanel(null);
        // panel.setEnabled(false);
        
        JFrame frame = new JFrame("FlagValuePanel Demo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(panel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
	
	/**
	 * Creates and return list of <code>DataTrimmingOperations</code> according to 
	 * given string.
	 * 
	 * @param discrete
	 * @return list of data trimming operations
	 */
	public static List<DataTrimmingOperation> parseDiscreteString(String discrete, String replacement, int columnIndex){
		List<DataTrimmingOperation> operations = new ArrayList<DataTrimmingOperation>();
		String[] splitted = discrete.split(";");
		for(String oldString : splitted){
			DataTrimmingOperation operation;
			try{
				Double stringAsNumber = Double.parseDouble(oldString);
				operation = new NumberToStringReplace(stringAsNumber, replacement, columnIndex);
			} catch (NumberFormatException nfe) {
				operation = new NormalStringReplace(oldString, replacement, true, columnIndex);
			}
			
			operations.add(operation);
		}
		return operations;
	}
}
