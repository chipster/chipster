
package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JToggleButton;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.ChipVisualisation;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.databeans.DataBean;

/**
 * Class for 3d scatterplot implementing functionality required in interface Visualisation. The
 * side panel is done and handled here, but the actual visualisation is drawn in CoordinateArea.
 * 
 * @author Petri Klemelä
 *
 */
public class Scatterplot3D extends ChipVisualisation implements ActionListener, KeyListener{
	
	public Scatterplot3D(VisualisationFrame frame) {
		super(frame);
	}

	private JPanel paramPanel;
    private JPanel settingsPanel;
    private AnnotateListPanel list;
    
    private JToggleButton rotateTool = new JToggleButton(VisualConstants.ROTATE_IMAGE);
    private JToggleButton handTool = new JToggleButton(VisualConstants.HAND_ICON);
    private JToggleButton selectTool = new JToggleButton(VisualConstants.ARROW_ICON);
    
    private JCheckBox autoCheckBox;
    
    public enum Tool {ROTATE, MOVE, SELECT};
    private Tool tool = Tool.ROTATE;
    
    private JButton toXY = new JButton(VisualConstants.XY_PLANE);
    private JButton toXZ = new JButton(VisualConstants.XZ_PLANE);
    private JButton toYZ = new JButton(VisualConstants.YZ_PLANE);
    
    private JComboBox xBox;
    private JComboBox yBox;
    private JComboBox zBox;
    private JComboBox colorBox;
    
    private JButton useButton;
    
    private final ClientApplication application = Session.getSession().getApplication();       
    
    private CoordinateArea coordinateArea;
    
    private static final Cursor ROTATE_CURSOR = Toolkit.getDefaultToolkit().
		createCustomCursor(VisualConstants.ROTATE_CURSOR_IMAGE.getImage(), 
				new Point(16,16), "Rotate");
    
    private static final Cursor ROTATE_AND_ZOOM_CURSOR = Toolkit.getDefaultToolkit().
	createCustomCursor(VisualConstants.ROTATE_AND_ZOOM_CURSOR_IMAGE.getImage(), 
			new Point(16,16), "Rotate");    
    
    private DataModel dataModel = new DataModel();
    
    private DataBean data;

	@Override
	public JPanel getParameterPanel() {
		//Recreate whole panel to make sure everything is updated if the dataset is changed
		if(paramPanel == null || data != 
			application.getSelectionManager().getSelectedDataBean()){
			
			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

			rotateTool.addActionListener(this);
			rotateTool.setToolTipText("Rotate by dragging, zoom with scroll wheel or with " +
					"Shift-button");
			handTool.addActionListener(this);
			handTool.setToolTipText("Move viewable part of the model by dragging it");
			selectTool.addActionListener(this);
			selectTool.setToolTipText("Selected datapoints are listed " +
					"in tab \'Selected\' below");

			rotateTool.addKeyListener(this);

			rotateTool.setSelected(true);

			toXY.addActionListener(this);
			toXY.setToolTipText("Show XY-plane");
			toXZ.addActionListener(this);		
			toXZ.setToolTipText("Show XZ-plane");
			toYZ.addActionListener(this);		
			toYZ.setToolTipText("Show YZ-plane");
			
			autoCheckBox = new JCheckBox("Automated rotation");
			autoCheckBox.addActionListener(this);

			JPanel settings = this.createSettingsPanel();
			list = new AnnotateListPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);
			tabPane.addTab("Selected", list);

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.insets.set(2, 10, 2, 10);
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.weighty = 0;
			c.weightx = 1.0;
			paramPanel.add(rotateTool, c);
			c.gridx++;
			paramPanel.add(handTool, c);
			c.gridx++;
			paramPanel.add(selectTool, c);
			c.gridy++;
			c.gridx = 0;
			paramPanel.add(toXY, c);
			c.gridx++;
			paramPanel.add(toXZ, c);
			c.gridx++;
			paramPanel.add(toYZ, c);

			c.gridx = 0;
			c.gridy++;
			c.gridwidth = 3;
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			paramPanel.add(autoCheckBox,c);
			c.gridy++;
			c.insets.set(5, 0, 0, 0);
			paramPanel.add(tabPane, c);
		}
		
		setToolsEnabled(false);
		autoCheckBox.setSelected(false);
		refreshAxisBoxes(application.getSelectionManager().getSelectedDataBean());
		
		return paramPanel;
	}
	
	public JPanel createSettingsPanel(){

		settingsPanel =  new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		xBox = new JComboBox();
		yBox = new JComboBox();
		zBox = new JComboBox();
		colorBox = new JComboBox();

		useButton = new JButton("Draw");
		useButton.addActionListener(this);

		JLabel xLabel = new JLabel("X-axis: ");
		xLabel.setForeground(Color.RED.darker().darker());
		JLabel yLabel = new JLabel("Y-axis: ");
		yLabel.setForeground(Color.GREEN.darker().darker());
		JLabel zLabel = new JLabel("Z-axis: ");
		zLabel.setForeground(Color.BLUE.darker().darker());
		JLabel colorLabel = new JLabel("Color: ");

		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.gridx = 0;
		c.insets.set(5, 10, 5, 10);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;
		c.gridx = 0;
		c.gridwidth = 3;
		settingsPanel.add(xLabel,c);		
		c.gridy++;
		settingsPanel.add(xBox,c);     
		c.gridy++;		
		settingsPanel.add(yLabel,c);
		c.gridy++;
		settingsPanel.add(yBox,c);
		c.gridy++;		
		settingsPanel.add(zLabel,c);
		c.gridy++;
		settingsPanel.add(zBox,c);
		c.gridy++;		
		settingsPanel.add(colorLabel,c);
		c.gridy++;
		settingsPanel.add(colorBox,c);
		c.gridy++;
		settingsPanel.add(useButton,c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;				
		settingsPanel.add(new JPanel(),c);
		
		return settingsPanel;
	}
	
	protected AnnotateListPanel getAnnotateList() {
		return list;
	}

	private void setToolsEnabled(boolean enabled) {
		rotateTool.setEnabled(enabled);
		handTool.setEnabled(enabled);
		selectTool.setEnabled(enabled);
		toXY.setEnabled(enabled);
		toXZ.setEnabled(enabled);
		toYZ.setEnabled(enabled);
	}

	private void refreshAxisBoxes(DataBean data) {
		if (paramPanel == null) {
			throw new IllegalStateException("must call getParameterPanel first");
		}
		this.updateCombo(xBox, data);
		this.updateCombo(yBox, data);
		this.updateCombo(zBox, data);
		this.updateCombo(colorBox, data);
	}
	
	private void updateCombo(JComboBox box, DataBean data){
		Visualisation.fillCompoBox(box, this.getVariablesFor(data));
	}
	
	public Tool getTool(){
		return tool;
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel
	 * to listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		
		
		
		if ( source == useButton ) {
			
			List<Variable> vars = new ArrayList<Variable>();
			vars.add((Variable)xBox.getSelectedItem());
			vars.add((Variable)yBox.getSelectedItem());
			vars.add((Variable)zBox.getSelectedItem());
			vars.add((Variable)colorBox.getSelectedItem());
			
			application.setVisualisationMethod(new VisualisationMethodChangedEvent(this,
					VisualisationMethod.SCATTERPLOT3D, vars, 
					getFrame().getDatas(), getFrame().getType(), getFrame()));
			
			coordinateArea.setPaintMode(CoordinateArea.PaintMode.PIXEL);			
			
		} else if (source == toXY){
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(Math.PI,0,0,1000,15);			
		} else if (source == toXZ){
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(-Math.PI/2+2*Math.PI,0,0,1000,15);
		} else if (source == toYZ){
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(Math.PI, -Math.PI / 2,0,1000,15);
		} 
		
		if(source == rotateTool || source == handTool || source == selectTool){
			this.selectTool((JToggleButton)source);			
		} 
		
		if(source == autoCheckBox){
			if(autoCheckBox.isSelected()){
				coordinateArea.movement.startAutomatedRotation();
			} else {
				stopAutoRotation();
			}
		}
	}
	
	public void stopAutoRotation(){
		autoCheckBox.setSelected(false);
		coordinateArea.movement.clearTasks();
	}
	
	private void selectTool(JToggleButton button){
		rotateTool.setSelected(button == rotateTool);
		handTool.setSelected(button == handTool);
		selectTool.setSelected(button == selectTool);
		
		if(button == rotateTool){
			coordinateArea.setCursor(ROTATE_CURSOR);
			tool = Tool.ROTATE;
		} else if(button == handTool){
			coordinateArea.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
			tool = Tool.MOVE;
		} else if(button == selectTool){
			coordinateArea.setCursor(Cursor.getDefaultCursor());
			tool = Tool.SELECT;
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		
		this.data = data;
		
		refreshAxisBoxes(data);	
		
		List<Variable> vars = getFrame().getVariables();
		if(vars == null || vars.size() < 4){
			if(xBox.getItemCount() >= 1){
				xBox.setSelectedIndex(0);
			}
			if(yBox.getItemCount() >= 2){
				yBox.setSelectedIndex(1);
			}
			if(zBox.getItemCount() >= 3){
				zBox.setSelectedIndex(2);
			}
			if(colorBox.getItemCount() >= 4){
				colorBox.setSelectedIndex(3);
			}						
		} else {		
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
			zBox.setSelectedItem(vars.get(2));
			colorBox.setSelectedItem(vars.get(3));
		}


		List<Variable> variables = new LinkedList<Variable>();
		variables.add((Variable)xBox.getSelectedItem()); 
		variables.add((Variable)yBox.getSelectedItem());
		variables.add((Variable)zBox.getSelectedItem());
		variables.add((Variable)colorBox.getSelectedItem());

		if (variables.size() >= 4 && 
				variables.get(0) != null && variables.get(1) != null &&
				variables.get(2) != null && variables.get(3) != null) {

			Iterable<String> identifier = data.queryFeatures("/column/ ").asStrings();
			Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
			Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
			Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();
			Iterable<Float> cValues = data.queryFeatures(variables.get(3).getExpression()).asFloats();
			
			dataModel.setData(identifier, xValues, yValues, zValues, cValues);

	    	JPanel panel = new JPanel();
	        panel.setLayout(new BorderLayout());
	        
	        coordinateArea = new CoordinateArea(this);
	        coordinateArea.addKeyListener(this);
	        coordinateArea.requestFocus();
	        panel.add(coordinateArea, BorderLayout.CENTER);
	        panel.add(new ColorScalePanel(dataModel), BorderLayout.WEST);
	        panel.setCursor(ROTATE_CURSOR);
	        this.setToolsEnabled(true);
	        
			return panel;
		}
		return this.getDefaultVisualisation();
	}
	
	public DataModel getDataModel() {
		return dataModel;
	}

	public void keyPressed(KeyEvent e) {					
		if(this.tool == Tool.ROTATE){
			if(e.isShiftDown()){
				coordinateArea.setCursor(ROTATE_AND_ZOOM_CURSOR);
			} else {
				coordinateArea.setCursor(ROTATE_CURSOR);
			}
		}
	}

	public void keyReleased(KeyEvent e) {
		if(this.tool == Tool.ROTATE){			
			coordinateArea.setCursor(ROTATE_CURSOR);			
		}
	}

	public void keyTyped(KeyEvent e) {		
		// ignore
	}
}

