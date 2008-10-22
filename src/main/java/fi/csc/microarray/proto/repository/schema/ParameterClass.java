package fi.csc.microarray.proto.repository.schema;

import java.util.ArrayList;
import java.util.List;


public class ParameterClass extends ParameterItem {
	
	private List<ParameterClass> subclasses;
	private List<ParameterInstance> instances;
	private int grade;
	
	public ParameterClass(String name) {
		super(name);
		this.subclasses = null;
		this.instances = null;
		this.grade = 0;
	}
	
	public void addSubclass(ParameterClass subclass) {
		if (subclass != null) {
			if (subclasses == null) {
				subclasses = new ArrayList<ParameterClass>();
			}
			subclasses.add(subclass);
			subclass.setParent(this);
		}
	}
	
	public ParameterClass[] getSubclasses() {
		if (subclasses != null) {
			return (ParameterClass[]) subclasses.toArray(new ParameterClass[] { });
		} else {
			return null;
		}
	}
	
	public int getSubclassCount() {
		if (subclasses != null) {
			return subclasses.size();
		} else {
			return 0;
		}
	}
	
	public boolean hasSubclasses() {
		return (subclasses != null && subclasses.size() > 0);
	}
	
	public void addInstance(ParameterInstance instance) {
		if (instance != null) {
			if (instances == null) {
				instances = new ArrayList<ParameterInstance>();
			}
			instances.add(instance);
			instance.setParent(this);
		}
	}
	
	public ParameterInstance[] getInstances() {
		if (instances != null) {
			return (ParameterInstance[]) instances.toArray(new ParameterInstance[] { });
		} else {
			return null;
		}
	}
	
	public int getInstanceCount() {
		if (instances != null) {
			return instances.size();
		} else {
			return 0;
		}
	}
	
	public boolean hasInstances() {
		return (instances != null && instances.size() > 0);
	}
	
	public int getGrade() {
		return grade;
	}
	
	public int setGrade(int newGrade) {
		this.grade = newGrade;
		int biggestGrade = this.grade;
		if (subclasses != null) {
			for (ParameterClass subclass : subclasses) {
				int biggestGradeWithin = subclass.setGrade(newGrade+1);
				if (biggestGradeWithin > biggestGrade) {
					biggestGrade = biggestGradeWithin;
				}
			}
		}
		return biggestGrade;
	}
	
	/**
	 * Returns null if feature is not supported or values are not available. 
	 */
	public Iterable<String> getPossibleValues() {
		return null;
	}
}