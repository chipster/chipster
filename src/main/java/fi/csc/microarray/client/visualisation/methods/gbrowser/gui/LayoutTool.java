package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

public class LayoutTool {
	
	public static enum LayoutMode { FIXED, FILL, FULL }

	public static void doLayout(LayoutContainer container, int layoutHeight) {
		int componentFillHeight = 0;
		
		if (inferLayoutMode(container) != LayoutMode.FIXED) {
			componentFillHeight = getComponentFillHeight(container, layoutHeight);
		}
		
		for (LayoutComponent component : container.getLayoutComponents()) {
			int childHeight;
			if (component.getLayoutMode() == LayoutMode.FIXED) {
				childHeight = component.getHeight();
			} else {
				childHeight = component.getMinHeight() + componentFillHeight;
				component.setHeight(childHeight);				
			}
				
			if (component instanceof LayoutContainer) {
				LayoutContainer childContainer = (LayoutContainer) component;
				doLayout(childContainer, childHeight);
//				component.setHeight(calculateHeight(childContainer));
			}
		}
	}
//	
//	public static void doCanvasLayout(LayoutContainer container, int layoutHeight) {
//		int nonFixedComponentHeight = 0;
//		
//		if (!isFixedHeight(container)) {
//			nonFixedComponentHeight = getComponentFillHeight(container, layoutHeight);
//		}
//		
//		for (LayoutComponent component : container.getLayoutComponents()) {
//			int childHeight;
//			if (component.isFixedHeight()) {
//				childHeight = component.getHeight();
//			} else {
//				childHeight = nonFixedComponentHeight;
//				childHeight = Math.max(childHeight, component.getMinHeight());
//				component.setHeight(childHeight);				
//			}
//				
//			if (component instanceof LayoutContainer) {
//				LayoutContainer childContainer = (LayoutContainer) component;
//				doViewPortLayout(childContainer, childHeight);
//				component.setHeight(calculateHeight(childContainer));
//			}
//		}
//	}

	private static int getComponentFillHeight(LayoutContainer container, int layoutHeight) {

		int fillHeight = layoutHeight - getMinHeightSum(container);
		return fillHeight / getFillComponentCount(container);
	}

	private static int getFillComponentCount(LayoutContainer container) {
		int count = 0;

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible() && component.getLayoutMode() != LayoutMode.FIXED) {
				count++;
			}
		}
		return count;	
	}

	public static int getMinHeightSum(LayoutContainer container) {
		int height = 0;

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible()) {
				if (component.getLayoutMode() == LayoutMode.FIXED) {
					height += component.getHeight();
				} else {
					height += component.getMinHeight();
				}
			}
		}
		return height;
	}

	/**
	 * @param container
	 * @return true if all components have fixed height, otherwise false
	 */
	public static LayoutMode inferLayoutMode(LayoutContainer container) {

		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible() && component.getLayoutMode() != LayoutMode.FIXED) {
				return LayoutMode.FILL;
			}
		}
		return LayoutMode.FIXED;
	}

	public static int getHeight(LayoutContainer container, int layoutHeight) {
		if (inferLayoutMode(container) == LayoutMode.FIXED) {
			return getMinHeightSum(container);
		} else {
			return layoutHeight;
		}
	}
//
//	public static int calculateHeight(LayoutContainer container) {
//		int height = 0;
//		for (LayoutComponent component : container.getLayoutComponents()) {
//			if (component.isVisible()) {
//				height += component.getHeight();
//			}
//		}
//		return height;
//	}

	public static int getFullHeight(LayoutContainer container) {
		int height = 0;
		for (LayoutComponent component : container.getLayoutComponents()) {
			if (component.isVisible()) {
				height += Math.max(component.getFullHeight(), component.getHeight());
			}
		}
		return height;		
	}
}
