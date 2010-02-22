
package fi.csc.microarray.client.visualisation.methods.gbrowser;

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
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JToggleButton;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.ChipVisualisation;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.fs.FSDataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;

/**
 * Class for 3d scatterplot implementing functionality required in interface Visualisation. The
 * side panel is done and handled here, but the actual visualisation is drawn in CoordinateArea.
 * 
 * @author Petri KlemelÃ¯Â¿Â½
 *
 */
public class GenomeBrowser extends Visualisation implements ActionListener, KeyListener{
	
	public GenomeBrowser(VisualisationFrame frame) {
		super(frame);
	}

	protected JPanel paramPanel;
    private JPanel settingsPanel;
    
    private JButton useButton;
    
    protected final ClientApplication application = Session.getSession().getApplication();       
    
    protected GenomePlot plot;

    protected DataBean data;
	private FileBrokerClient fileBroker;

	@Override
	public JPanel getParameterPanel() {

		if(paramPanel == null || data != 
			application.getSelectionManager().getSelectedDataBean()){
			
			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			
			
			paramPanel.add(new JLabel("Chromosome"), c);
			
			JPanel settings = this.createSettingsPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			c.weightx = 1.0;
			c.insets.set(5, 0, 0, 0);
			
			paramPanel.add(tabPane, c);
			
			
			
		}				
		
		return paramPanel;
	}
	
	public void createAnnotationComponents(JPanel panel, GridBagConstraints c) {
		
		MessagingEndpoint endpoint;
		try {
			endpoint = new MessagingEndpoint(new NodeBase() {
				public String getName() {
					return "gbrowserAnnotations";
				}
			});
			
			this.fileBroker = new FileBrokerClient(endpoint.createTopic(Topics.Name.URL_TOPIC, AccessMode.READ));
			
			fileBroker.getFile(new URL("http://chipster-devel.csc.fi:8080/public/ensembl/contents.txt"));
			
		} catch (MicroarrayException e) {
			application.reportException(e);
		} catch (JMSException e) {
			application.reportException(e);
		} catch (MalformedURLException e) {
			application.reportException(e);
		} catch (IOException e) {
			application.reportException(e);
		}
	}
	
	public JPanel createSettingsPanel(){

		settingsPanel =  new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		useButton = new JButton("Draw");
		useButton.addActionListener(this);
		
		JComboBox chrBox = new JComboBox();
		
		//FIXME These should be read from user data file
		for( int i = 1; i <= 22; i++) {
			chrBox.addItem(i);
		}
		
		chrBox.addActionListener(this);
		
		JTextArea megaLocation = new JTextArea(1, 3);
		JTextArea kiloLocation = new JTextArea(1, 3);
		JTextArea unitLocation = new JTextArea(1, 3);

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
				
		settingsPanel.add(new JLabel("Chromosome"));		
		c.gridy++;
		settingsPanel.add(chrBox, c);     
		c.gridy++;		
		settingsPanel.add(new JLabel("Location"),c);
		c.gridy++;
		settingsPanel.add(megaLocation, c);
		c.gridx++;
		settingsPanel.add(new JLabel("M"), c);
		c.gridx++;		
		settingsPanel.add(kiloLocation, c);
		c.gridx++;
		settingsPanel.add(new JLabel("k"), c);
		c.gridx++;
		settingsPanel.add(unitLocation, c);
		c.gridx = 0;
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
	
	protected JComponent getColorLabel() {
		return new JLabel("Color: ");
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

	protected void refreshAxisBoxes(DataBean data) {
		if (paramPanel == null) {
			throw new IllegalStateException("must call getParameterPanel first");
		}
		this.updateCombo(xBox, data);
		this.updateCombo(yBox, data);
		this.updateCombo(zBox, data);
		Visualisation.fillCompoBox(colorBox, this.getVariablesMore(data));
	}
	
	protected void updateCombo(JComboBox box, DataBean data){
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
			
			useButtonPressed();
			
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
	
	protected void useButtonPressed() {
		List<Variable> vars = new ArrayList<Variable>();
		vars.add((Variable)xBox.getSelectedItem());
		vars.add((Variable)yBox.getSelectedItem());
		vars.add((Variable)zBox.getSelectedItem());
		vars.add((Variable)colorBox.getSelectedItem());
		
		application.setVisualisationMethod(new VisualisationMethodChangedEvent(this,
				VisualisationMethod.SCATTERPLOT3D, vars, 
				getFrame().getDatas(), getFrame().getType(), getFrame()));
		
		coordinateArea.setPaintMode(CoordinateArea.PaintMode.PIXEL);
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

			retrieveData(variables);

	    	JPanel panel = new JPanel();
	        panel.setLayout(new BorderLayout());
	        
	        coordinateArea = new CoordinateArea(this);
	        coordinateArea.addKeyListener(this);
	        coordinateArea.requestFocus();
	        
	        JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	        split.setDividerSize(3);
	        split.setOpaque(true);
	        split.setRightComponent(coordinateArea);
	        split.setLeftComponent(getColorScalePanel());
	        split.setContinuousLayout(true);
	        
	        panel.add(split, BorderLayout.CENTER);
	        
//	        panel.add(coordinateArea, BorderLayout.CENTER);
//	        panel.add(getColorScalePanel(), BorderLayout.WEST);
	        
	        coordinateArea.setCursor(ROTATE_CURSOR);
	        this.setToolsEnabled(true);
	        
			return panel;
		}
		return this.getDefaultVisualisation();
	}
	
	protected JComponent getColorScalePanel(){
		return new ColorScalePanel(dataModel);
	}
	
	
	protected void retrieveData(List<Variable> variables) throws MicroarrayException {
		Iterable<String> identifier = data.queryFeatures("/identifier").asStrings();
		Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
		Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();
		Iterable<Float> cValues = data.queryFeatures(variables.get(3).getExpression()).asFloats();
		
		dataModel.setData(identifier, xValues, yValues, zValues, cValues);
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
	
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return super.canVisualise(bean) && 
			!VisualisationMethod.SCATTERPLOT3DPCA.getHeadlessVisualiser().canVisualise(bean);
		
	}
}

