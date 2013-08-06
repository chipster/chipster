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
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JToggleButton;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.methods.ChipVisualisation;
import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea.PaintMode;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.ScaleUtil;

/**
 * Class for 3d scatterplot implementing functionality required in interface Visualisation. The side panel is done and handled here, but the
 * actual visualisation is drawn in CoordinateArea.
 * 
 * @author Petri Klemelä
 * 
 */
public class Scatterplot3D extends ChipVisualisation implements ActionListener, KeyListener {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	protected JPanel paramPanel;
	private JPanel settingsPanel;
	private AnnotateListPanel list;

	private JToggleButton rotateTool = new JToggleButton(VisualConstants.ROTATE_IMAGE);
	private JToggleButton handTool = new JToggleButton(VisualConstants.HAND_ICON);
	private JToggleButton selectTool = new JToggleButton(VisualConstants.ARROW_ICON);

	private JCheckBox autoCheckBox;

	public enum Tool {
		ROTATE, MOVE, SELECT
	};

	private Tool tool = Tool.ROTATE;

	private JButton toXY = new JButton(VisualConstants.XY_PLANE);
	private JButton toXZ = new JButton(VisualConstants.XZ_PLANE);
	private JButton toYZ = new JButton(VisualConstants.YZ_PLANE);

	protected JComboBox xBox;
	protected JComboBox yBox;
	protected JComboBox zBox;
	protected JComboBox colorBox;
	
	private JCheckBox blackBackgroundCheckBox;
	private JLabel paintModeLabel;
	private JComboBox<PaintMode> paintModeBox;

	private JButton useButton;

	protected final ClientApplication application = Session.getSession().getApplication();

	protected CoordinateArea coordinateArea;

	private static final Cursor ROTATE_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(VisualConstants.ROTATE_CURSOR_IMAGE.getImage(), new Point(16, 16), "Rotate");

	private static final Cursor ROTATE_AND_ZOOM_CURSOR = Toolkit.getDefaultToolkit().createCustomCursor(VisualConstants.ROTATE_AND_ZOOM_CURSOR_IMAGE.getImage(), new Point(16, 16), "Rotate");
	private static final int DEFAULT_TO_DOT_PAINT_MODE = 20000;

	protected DataBean data;
	protected ColorScalePanel scalePanel;
	protected DataModel dataModel;

	@Override
	public JPanel getParameterPanel() {
		// Recreate whole panel to make sure everything is updated if the dataset is changed
		if (paramPanel == null || data != application.getSelectionManager().getSelectedDataBean()) {

			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

			rotateTool.addActionListener(this);
			rotateTool.setToolTipText("Rotate by dragging, zoom with scroll wheel or with " + "Shift-button");
			handTool.addActionListener(this);
			handTool.setToolTipText("Move viewable part of the model by dragging it");
			selectTool.addActionListener(this);
			selectTool.setToolTipText("Selected datapoints are listed " + "in tab \'Selected\' below");

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
			
			blackBackgroundCheckBox = new JCheckBox("Black background");
			blackBackgroundCheckBox.setSelected(false);
			blackBackgroundCheckBox.addActionListener(this);
			
			paintModeLabel = new JLabel("Symbol: ");
			paintModeBox = new JComboBox<PaintMode>(PaintMode.values());
			paintModeBox.addActionListener(this);

			JPanel settings = this.createSettingsPanel();
			list = createListPanel();

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
			c.gridy++;
			c.gridx = 0;
			c.gridwidth = 3;
			paramPanel.add(autoCheckBox, c);
			c.gridy++;
			paramPanel.add(blackBackgroundCheckBox, c);
			c.gridy++;
			paramPanel.add(paintModeLabel, c);
			c.gridy++;
			paramPanel.add(paintModeBox, c);
			c.gridy++;
			paramPanel.add(new JLabel(" "), c);
			c.gridy++;
			
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			c.insets.set(5, 0, 0, 0);
			paramPanel.add(tabPane, c);
		}

		setToolsEnabled(false);
		autoCheckBox.setSelected(false);
		refreshAxisBoxes(application.getSelectionManager().getSelectedDataBean());

		return paramPanel;
	}

	/**
	 * Hook for PCA to change the word "genes"
	 */
	public AnnotateListPanel createListPanel() {
		return new AnnotateListPanel();
	}

	public JPanel createSettingsPanel() {

		settingsPanel = new JPanel();
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
		JComponent colorLabel = getColorLabel();

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
		settingsPanel.add(xLabel, c);
		c.gridy++;
		settingsPanel.add(xBox, c);
		c.gridy++;
		settingsPanel.add(yLabel, c);
		c.gridy++;
		settingsPanel.add(yBox, c);
		c.gridy++;
		settingsPanel.add(zLabel, c);
		c.gridy++;
		settingsPanel.add(zBox, c);
		c.gridy++;
		settingsPanel.add(colorLabel, c);
		c.gridy++;
		settingsPanel.add(colorBox, c);		
		c.gridy++;
		settingsPanel.add(useButton, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		settingsPanel.add(new JPanel(), c);

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
		Visualisation.fillComboBox(colorBox, this.getVariablesMore(data));
	}

	protected void updateCombo(JComboBox box, DataBean data) {
		Visualisation.fillComboBox(box, this.getVariablesFor(data));
	}

	public Tool getTool() {
		return tool;
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == useButton) {

			useButtonPressed();

		} else if (source == toXY) {
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(Math.PI, 0, 0, 1000, 15);
		} else if (source == toXZ) {
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(-Math.PI / 2 + 2 * Math.PI, 0, 0, 1000, 15);
		} else if (source == toYZ) {
			stopAutoRotation();
			coordinateArea.movement.addRotationTask(Math.PI, -Math.PI / 2, 0, 1000, 15);
		}

		if (source == rotateTool || source == handTool || source == selectTool) {
			this.selectTool((JToggleButton) source);
		}

		if (source == autoCheckBox) {
			if (autoCheckBox.isSelected()) {
				coordinateArea.movement.startAutomatedRotation();
			} else {
				stopAutoRotation();
			}
		}
		
		if (source == blackBackgroundCheckBox) {
			
			updateBackgroundColor();
		}
		
		if (source == paintModeBox) {
			
			coordinateArea.setPaintMode((PaintMode) paintModeBox.getSelectedItem());
		}
	}
	
	private void updateBackgroundColor() {
		Color bg;
		Color fg;
		
		if (blackBackgroundCheckBox.isSelected()) {
			bg = Color.black;
			fg = Color.white;
		} else {
			bg = Color.white;
			fg = Color.black;
		}			
		coordinateArea.setBackground(bg);
		coordinateArea.setForeground(fg);
		scalePanel.setBackground(bg);
		scalePanel.setForeground(fg);
	}

	protected void useButtonPressed() {
		List<Variable> vars = new ArrayList<Variable>();
		vars.add((Variable) xBox.getSelectedItem());
		vars.add((Variable) yBox.getSelectedItem());
		vars.add((Variable) zBox.getSelectedItem());
		vars.add((Variable) colorBox.getSelectedItem());

		application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, MicroarrayModule.VisualisationMethods.SCATTERPLOT3D, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));
	}

	public void stopAutoRotation() {
		autoCheckBox.setSelected(false);
		coordinateArea.movement.clearTasks();
	}

	private void selectTool(JToggleButton button) {
		rotateTool.setSelected(button == rotateTool);
		handTool.setSelected(button == handTool);
		selectTool.setSelected(button == selectTool);

		if (button == rotateTool) {
			coordinateArea.setCursor(ROTATE_CURSOR);
			tool = Tool.ROTATE;
		} else if (button == handTool) {
			coordinateArea.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
			tool = Tool.MOVE;
		} else if (button == selectTool) {
			coordinateArea.setCursor(Cursor.getDefaultCursor());
			tool = Tool.SELECT;
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		refreshAxisBoxes(data);

		List<Variable> vars = getFrame().getVariables();
		if (vars == null || vars.size() < 4) {
			if (xBox.getItemCount() >= 1) {
				xBox.setSelectedIndex(0);
			}
			if (yBox.getItemCount() >= 2) {
				yBox.setSelectedIndex(1);
			}
			if (zBox.getItemCount() >= 3) {
				zBox.setSelectedIndex(2);
			}

			if (colorBox.getItemCount() >= 4) {
				colorBox.setSelectedIndex(3);
			}
		} else {
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
			zBox.setSelectedItem(vars.get(2));
			colorBox.setSelectedItem(vars.get(3));
		}

		List<Variable> variables = new LinkedList<Variable>();
		variables.add((Variable) xBox.getSelectedItem());
		variables.add((Variable) yBox.getSelectedItem());
		variables.add((Variable) zBox.getSelectedItem());
		variables.add((Variable) colorBox.getSelectedItem());

		if (variables.size() >= 4 && variables.get(0) != null && variables.get(1) != null && variables.get(2) != null && variables.get(3) != null) {

			retrieveData(variables);
			

			JPanel panel = new JPanel();
			panel.setLayout(new BorderLayout());

			coordinateArea = new CoordinateArea(this);
			coordinateArea.addKeyListener(this);
			coordinateArea.requestFocus();
			
			if (getDataModel().getDataArray().length > DEFAULT_TO_DOT_PAINT_MODE) {
				paintModeBox.setSelectedItem(PaintMode.RECT);
			}

			JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
			split.setDividerSize(3);
			split.setOpaque(true);
			split.setRightComponent(coordinateArea);
			scalePanel = null;
			split.setLeftComponent(getColorScalePanel());
			split.setContinuousLayout(true);
			split.setDividerLocation(150);

			panel.add(split, BorderLayout.CENTER);

			coordinateArea.setCursor(ROTATE_CURSOR);
			this.setToolsEnabled(true);
			
			updateBackgroundColor();

			return panel;
		}
		return this.getDefaultVisualisation();
	}

	protected JComponent getColorScalePanel() {
		
		if (scalePanel == null) {
			List<String> colorGroupNames = new LinkedList<>();
			List<Float> colorScaleValues = new LinkedList<>();

			float[] scale = getDataModel().getColorScaleValues();
			for (int i = 0; i < scale.length - 1; i++) {
				String floor = ScaleUtil.format(scale[i]);
				String ceiling = ScaleUtil.format(scale[i+1]);
				colorGroupNames.add("" + floor + " - " + ceiling);
				colorScaleValues.add(scale[i]);
			}

			scalePanel = new ColorScalePanel(getDataModel(), colorGroupNames, colorScaleValues);
		}
		return scalePanel;
	}

	protected void retrieveData(List<Variable> variables) throws MicroarrayException {
		Iterable<String> identifier = data.queryFeatures("/identifier").asStrings();
		Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
		Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();
		Iterable<Float> cValues = data.queryFeatures(variables.get(3).getExpression()).asFloats();		

		getDataModel().setData(identifier, xValues, yValues, zValues, cValues);	
	}

	public DataModel getDataModel() {
		
		if (dataModel == null) {
			dataModel = new DataModel();
		}
		return dataModel;
	}

	public void keyPressed(KeyEvent e) {
		if (this.tool == Tool.ROTATE) {
			if (e.isShiftDown()) {
				coordinateArea.setCursor(ROTATE_AND_ZOOM_CURSOR);
			} else {
				coordinateArea.setCursor(ROTATE_CURSOR);
			}
		}
	}

	public void keyReleased(KeyEvent e) {
		if (this.tool == Tool.ROTATE) {
			coordinateArea.setCursor(ROTATE_CURSOR);
		}
	}

	public void keyTyped(KeyEvent e) {
		// ignore
	}
}
