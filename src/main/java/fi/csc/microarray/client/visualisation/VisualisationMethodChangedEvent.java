package fi.csc.microarray.client.visualisation;

import java.beans.PropertyChangeEvent;
import java.util.List;

import fi.csc.microarray.client.visualisation.VisualisationFactory.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.databeans.DataBean;


public class VisualisationMethodChangedEvent extends PropertyChangeEvent {
	
	private List<DataBean> datas;
	private List<Variable> variables;
	private VisualisationFrame targetFrameInstance;
	
	public List<DataBean> getDatas() {
		return datas;
	}

	public FrameType getTarget() {
		return target;
	}
	
	public List<Variable> getVariables(){
		return variables;
	}
	
	/**
	 * FrameInstance can be used to forward the visualisation events to the specific external window, 
	 * because they all share the same target type
	 * 
	 * @return
	 */
	public VisualisationFrame getTargetFrameInstance(){
		return targetFrameInstance;
	}

	private FrameType target;

	public VisualisationMethodChangedEvent(Object source, VisualisationMethod newMethod, 
			List<Variable> variables, List<DataBean> datas, FrameType target ) {
		
		super(source, null, null, newMethod);
		this.datas = datas;
		this.target = target;
		this.variables = variables;
	}
	
	/**
	 * FrameInstance can be used to forward the visualisation events to the specific external window, 
	 * because they all share the same target type.
	 * 
	 * @param source
	 * @param newMethod
	 * @param variables
	 * @param datas
	 * @param target
	 * @param frameInstance
	 */
	public VisualisationMethodChangedEvent(Object source, VisualisationMethod newMethod, 
			List<Variable> variables, List<DataBean> datas, FrameType target, 
			VisualisationFrame frameInstance ) {
		
		this(source, newMethod, variables, datas, target);
		this.targetFrameInstance = frameInstance;
	}

	public VisualisationMethod getNewMethod() {
		return (VisualisationMethod)super.getNewValue();
	}

	public void setTargetFrameInstance(ExternalVisualisationFrame window) {
		targetFrameInstance = window;
	}		
}
