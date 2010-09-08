package fi.csc.microarray.client.visualisation.methods;

import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.threed.ColorGroupsPanel;
import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class Scatterplot3DPCA extends Scatterplot3D{
	
	DataBean phenoBean;
	private List<String> colorGroupList;
	private JScrollPane legendScroller;

	public Scatterplot3DPCA(VisualisationFrame frame) {
		super(frame);
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
		
		if(phenoBean != null){
			phenoCols = Arrays.asList(
					VisualisationUtilities.getVariablesFilteredInclusive(phenoBean, "", false));
		} else {
			phenoCols = new LinkedList<Variable>();
			phenoCols.add(new Variable("No phenodata", ""));
		}
		
		Visualisation.fillCompoBox(colorBox, phenoCols.toArray(new Variable[0]));
	}
	
	@Override
	protected void retrieveData(List<Variable> variables) throws MicroarrayException {
		Iterable<String> identifier = data.queryFeatures("/identifier").asStrings();
		Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
		Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();
		
		String colorQuery = variables.get(3).getExpression();;
		
		Iterable<String> colorStrings;
		
		if(phenoBean != null) {
			colorStrings = phenoBean.queryFeatures(colorQuery).asStrings();
		} else {
			LinkedList<String> list = new LinkedList<String>();
					
			Iterator<Float> iter = xValues.iterator();
			while(iter.hasNext()) {
				iter.next();
				list.add("");
			}
			colorStrings = list;
		}
		
		TreeSet<String> groups = new TreeSet<String>();
		for(String str : colorStrings){
			groups.add(str);			
		}
		
		colorGroupList = Arrays.asList(groups.toArray(new String[0]));
		
		List<Float> cValues = new LinkedList<Float>();
		
		for(String str : colorStrings) {
			cValues.add((float)colorGroupList.indexOf(str));
		}
		
		dataModel.setData(identifier, xValues, yValues, zValues, cValues);
	}
	
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		
		boolean isTabular = BasicModule.VisualisationMethods.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
		boolean isChips = false;
		Parameter pcaOn = bean.getOperation().getParameter("do.pca.on");
		if (pcaOn != null){
			Object value = pcaOn.getValue();
			if( value != null){
				if(value.equals("chips")){
					isChips = true;
				}
			}
		}
		return bean.getOperation().getDefinition().getID().equals("PCA") &&
		isChips &&
		isTabular && 
		hasRows(bean) && 
		bean.queryFeatures("/column/chip.*").exists();		
		}
	
	protected void useButtonPressed() {
		List<Variable> vars = new ArrayList<Variable>();
		vars.add((Variable)xBox.getSelectedItem());
		vars.add((Variable)yBox.getSelectedItem());
		vars.add((Variable)zBox.getSelectedItem());
		vars.add((Variable)colorBox.getSelectedItem());
		
		application.setVisualisationMethod(new VisualisationMethodChangedEvent(this,
				MicroarrayModule.VisualisationMethods.SCATTERPLOT3DPCA, vars, 
				getFrame().getDatas(), getFrame().getType(), getFrame()));
		
		coordinateArea.setPaintMode(CoordinateArea.PaintMode.PIXEL);
	}
	
	@Override
	protected JComponent getColorLabel() {
		
		JTextArea label = new JTextArea("Color from\nphenodata:");
		label.setOpaque(false);
		//label.setText("Color from\nphenodata:");
		label.setEditable(false);
		
		return label;
	}
	
	protected JComponent getColorScalePanel(){
		ColorGroupsPanel legendPanel = new ColorGroupsPanel(dataModel, colorGroupList);
		//legendPanel.setPreferredSize(new Dimension(300,300));
		legendScroller = new JScrollPane(legendPanel);
		legendScroller.setPreferredSize(new Dimension(50,300));
		//legendScroller.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		return legendScroller;
	}
}
