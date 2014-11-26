package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Component;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;

public class LayoutUpdater<ChildType extends LayoutChild, ParentType extends LayoutParent > {

	/**
	 * Update layout components. It would be much easier to just remove all old
	 * components and ad the new components, but this breaks mouse dragged
	 * events on OS X. Instead, this updates the old components in place and
	 * removes and adds the components only when necessary.
	 * 
	 * @param children
	 * @param parent
	 */
	public void update(List<ChildType> children, ParentType parent) {
		
		for (ChildType child : children) {
			child.updateLayout();
		}

		// create a new list to make it editable
		Component[] array = parent.getLayoutContainer().getComponents();
		List<Component> oldComponents = new ArrayList<>(Arrays.asList(array));
		
		// collect a list of new components
		List<Component> newComponents = new ArrayList<>();
		for (ChildType child : children) {
			Component newComponent = child.getLayoutComponent();
			newComponents.add(newComponent);
		}
		
		// remove components that aren't needed anymore
		for (Component oldComponent : oldComponents) {
			if (!newComponents.contains(oldComponent)) {
				parent.getLayoutContainer().remove(oldComponent);
			}
		}
		
		// check if old and new components are the same
		boolean keepOld = true;
		for (int i = 0; i < children.size(); i++) {
			ChildType child = children.get(i);
			Component newComponent = child.getLayoutComponent();
			Component oldComponent = null;
			if (parent.getLayoutContainer().getComponentCount() > i) {
				oldComponent = parent.getLayoutContainer().getComponent(i);
			}
			if (newComponent != oldComponent) {
				keepOld = false;
				break;
			}
		}
		
		if (!keepOld) {
			parent.getLayoutContainer().removeAll();
		}
		
		// update or add new components
		for (int i = 0; i < children.size(); i++) {
			ChildType child = children.get(i);
			Component newComponent = child.getLayoutComponent();
		
			if (keepOld) {
				MigLayout migLayout = (MigLayout)parent.getLayoutContainer().getLayout();
				migLayout.setComponentConstraints(newComponent, getConstraint(child.getLayoutMode()));
			} else {
				parent.getLayoutContainer().add(newComponent, getConstraint(child.getLayoutMode()));
			}
		}
	}

	public String getConstraint(LayoutMode mode) {
		String constraints;
		if (LayoutMode.FILL == mode) {
			constraints = "grow";
		} else {
			constraints =  "growx";
		}
		return constraints;
	}
}
