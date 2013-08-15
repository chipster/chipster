package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.threed.ColorScalePanel;
import fi.csc.microarray.client.visualisation.methods.threed.DataModel;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class Scatterplot3DPCA extends Scatterplot3D {
	
	public class PCADataModel extends DataModel {
		public PCADataModel(Color[] colorScheme) {
			super(colorScheme);
		}

		@Override
		public Color getColorFor(float value) {
			Color c = Color.gray;		
			int index = (int)value;
			if (index >= 0 && index < colorGroupNames.size() && index < getColorScheme().length) {
				return getColorScheme()[index];
			}
			return c;
		}
	}

	DataBean phenoBean;
	private List<String> colorGroupNames;
	private JScrollPane legendScroller;
	private LinkedList<Float> colorGroupValues;	

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public AnnotateListPanel createListPanel() {
		return new AnnotateListPanel("Chips", false);
	}

	@Override
	protected void refreshAxisBoxes(DataBean data) {
		if (paramPanel == null) {
			throw new IllegalStateException("must call getParameterPanel first");
		}
		this.updateCombo(xBox, data);
		this.updateCombo(yBox, data);
		this.updateCombo(zBox, data);

		phenoBean = LinkUtils.retrieveInherited(data, Link.ANNOTATION);

		List<Variable> phenoCols;

		if (phenoBean != null) {
			phenoCols = Arrays.asList(VisualisationUtilities.getVariablesFilteredInclusive(phenoBean, "", false));
		} else {
			phenoCols = new LinkedList<Variable>();
			phenoCols.add(new Variable("No phenodata", ""));
		}

		Visualisation.fillComboBox(colorBox, phenoCols.toArray(new Variable[0]));
	}

	@Override
	protected void retrieveData(List<Variable> variables) throws MicroarrayException {
		Iterable<String> identifier = data.queryFeatures("/identifier").asStrings();
		Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
		Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();

		String colorQuery = variables.get(3).getExpression();

		Iterable<String> colorDataStrings;

		if (phenoBean != null) {
			colorDataStrings = phenoBean.queryFeatures(colorQuery).asStrings();
		} else {
			LinkedList<String> list = new LinkedList<String>();

			Iterator<Float> iter = xValues.iterator();
			while (iter.hasNext()) {
				iter.next();
				list.add("");
			}
			colorDataStrings = list;
		}

		//retains only unique names
		Set<String> nameSet = new HashSet<>(); 		
		for (String str : colorDataStrings) {
			nameSet.add(str);
		}
		colorGroupNames = new LinkedList<String>(nameSet);
		Collections.sort(colorGroupNames);
		
		//This list contains example value for each unique name. 
		//The color of the group in ScalePanel is queried with this value  
		colorGroupValues = new LinkedList<Float>();
		
		for (String str : colorGroupNames) {
			colorGroupValues.add((float) colorGroupNames.indexOf(str));
		}
		
		//color index for each data point
		List<Float> colorDataValues = new LinkedList<Float>();
		for (String str : colorDataStrings) {
			colorDataValues.add((float)colorGroupNames.indexOf(str));
		}		

		getDataModel().setData(identifier, xValues, yValues, zValues, colorDataValues);
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return super.canVisualise(bean) && bean.hasTypeTag(MicroarrayModule.TypeTags.EXPRESSION_PRIMARY_COMPONENTS_CHIPWISE);
	}

	protected void useButtonPressed() {
		List<Variable> vars = new ArrayList<Variable>();
		vars.add((Variable) xBox.getSelectedItem());
		vars.add((Variable) yBox.getSelectedItem());
		vars.add((Variable) zBox.getSelectedItem());
		vars.add((Variable) colorBox.getSelectedItem());

		application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, MicroarrayModule.VisualisationMethods.SCATTERPLOT3DPCA, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));				
	}

	@Override
	protected JComponent getColorLabel() {

		JTextArea label = new JTextArea("Color from\nphenodata:");
		label.setOpaque(false);
		label.setEditable(false);

		return label;
	}

	@Override
	protected JComponent getColorScalePanel() {
		if (scalePanel == null) {
			scalePanel = new ColorScalePanel(getDataModel(), colorGroupNames, colorGroupValues);
			legendScroller = new JScrollPane(scalePanel);
			legendScroller.setPreferredSize(new Dimension(50, 300));
		}
		return legendScroller;
	}
	
	@Override
	public DataModel getDataModel() {
		
		if (dataModel == null) {
			dataModel = new PCADataModel(DataModel.qualitatativeColorScheme);
		}
		return dataModel;
	}
}
