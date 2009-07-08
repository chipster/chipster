package fi.csc.microarray.client.visualisation.methods;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.threed.CoordinateArea;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;

public class Scatterplot3DPCA extends Scatterplot3D{
	
	DataBean phenoBean;

	public Scatterplot3DPCA(VisualisationFrame frame) {
		super(frame);
	}
	
	@Override
	public AnnotateListPanel createListPanel() {
		return new AnnotateListPanel("Chips");
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
		List<Variable> phenoCols = Arrays.asList(VisualisationUtilities.getVariablesFilteredInclusive(phenoBean, "", false));
		List<Variable> colsToRemove = new ArrayList<Variable>();
		
		for(Variable col : phenoCols){
			
			Iterator<String> values = null;
			try {
				values = phenoBean.queryFeatures(col.getExpression()).asStrings().iterator();
				
				while(values.hasNext()){
					try {
						Integer.parseInt(values.next());
					} catch (NumberFormatException e) {
						colsToRemove.add(col);
						break;
					}
					
				}
			} catch (MicroarrayException e1) {
				application.reportException(new MicroarrayException(
						"Unable to find phenodata, try normal 3D Scatterplot instead", e1));
			}				
		}
		//java.lang.UnsupportedOperationException:
		//phenoCols.removeAll(colsToRemove);
		
		List<Variable> numericCols = new ArrayList<Variable>();
		
		for(Variable col: phenoCols){
			if(!colsToRemove.contains(col)){
				numericCols.add(col);
			}				
		}
		
		Visualisation.fillCompoBox(colorBox, numericCols.toArray(new Variable[0])); 				
	}
	
	@Override
	protected void retrieveData(List<Variable> variables) throws MicroarrayException {
		Iterable<String> identifier = data.queryFeatures("/identifier").asStrings();
		Iterable<Float> xValues = data.queryFeatures(variables.get(0).getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(variables.get(1).getExpression()).asFloats();
		Iterable<Float> zValues = data.queryFeatures(variables.get(2).getExpression()).asFloats();
		Iterable<Float> cValues = phenoBean.queryFeatures(variables.get(3).getExpression()).asFloats();
		
		dataModel.setData(identifier, xValues, yValues, zValues, cValues);
	}
	
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		
		boolean isTabular = VisualisationMethod.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
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
		return bean.getOperation().getDefinition().getName().equals("PCA") &&
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
				VisualisationMethod.SCATTERPLOT3DPCA, vars, 
				getFrame().getDatas(), getFrame().getType(), getFrame()));
		
		coordinateArea.setPaintMode(CoordinateArea.PaintMode.PIXEL);
	}

}
