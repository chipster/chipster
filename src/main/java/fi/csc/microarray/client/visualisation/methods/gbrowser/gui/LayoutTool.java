package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

public class LayoutTool {

	public static void doLayout(LayoutContainer container, int layoutHeight) {
		int nonFixedComponentHeight = 0;
		
		if (!isFixedHeight(container)) {
			nonFixedComponentHeight = getNonFixedComponentHeight(container, layoutHeight);
		}
		
		for (LayoutComponent component : container.getLayoutComponents()) {
			int childHeight;
			if (component.isFixedHeight()) {
				childHeight = component.getHeight();
			} else {
				childHeight = nonFixedComponentHeight;
				childHeight = Math.max(childHeight, component.getMinHeight());
				component.setHeight(childHeight);				
			}
				
			if (component instanceof LayoutContainer) {
				LayoutContainer childContainer = (LayoutContainer) component;
				doLayout(childContainer, childHeight);
				component.setHeight(calculateHeight(childContainer));
			}
		}
	}

	private static int getNonFixedComponentHeight(LayoutContainer container, int layoutHeight) {

		int nonFixedHeightSum = layoutHeight - getFixedHeightSum(container);
		return nonFixedHeightSum / getNonFixedHeightComponentCount(container);
	}

	private static int getNonFixedHeightComponentCount(LayoutContainer container) {
		int count = 0;

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible() && !component.isFixedHeight()) {
				count++;
			}
		}
		return count;	
	}

	private static int getFixedHeightSum(LayoutContainer container) {
		int height = 0;

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible() && component.isFixedHeight()) {
				height += component.getHeight();
			}
		}
		return height;
	}

	/**
	 * @param container
	 * @return true if all components have fixed height, otherwise false
	 */
	public static boolean isFixedHeight(LayoutContainer container) {

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible() && !component.isFixedHeight()) {
				return false;
			}
		}
		return true;
	}

	public static int getHeight(LayoutContainer container, int layoutHeight) {
		if (isFixedHeight(container)) {
			return getFixedHeightSum(container);
		} else {
			return layoutHeight;
		}
	}

	public static int calculateHeight(LayoutContainer container) {
		int height = 0;
		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible()) {
				height += component.getHeight();
			}
		}
		return height;
	}

	public static int getCanvasHeight(LayoutContainer container) {
		int height = 0;
		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible()) {
				height += Math.max(component.getCanvasHeight(), component.getHeight());
			}
		}
		return height;		
	}
}
