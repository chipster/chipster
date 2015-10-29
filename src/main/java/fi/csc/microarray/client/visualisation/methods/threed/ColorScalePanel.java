package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.RenderingHints;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

/**
 * @author Petri Klemelä
 */
public class ColorScalePanel extends JPanel{
	
	public static class BallIcon extends ImageIcon {

		private Color color;

		public BallIcon(Color color) {
			this.color = color;
		}
		
		@Override
		public void paintIcon(Component component, Graphics g, int x, int y) {
			super.paintIcon(component, g, x, y);
			Graphics2D g2 = (Graphics2D)g;
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON);

			DataPoint.paintBall(2, 2, 18, 18, color, g2);
		}
		
	    public int getIconWidth() {
	        return 24;
	    }
	    
	    public int getIconHeight() {
	        return 24;
	    }
	}

	DataModel dataModel;
	List<Float> colorGroupValues;

	public ColorScalePanel(DataModel dataModel, List<String> colorGroupNames, List<Float> colorGroupValues){
		this.dataModel = dataModel;
		this.colorGroupValues = colorGroupValues;
						
		this.setLayout(new GridBagLayout());

		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;

		for(String groupStr : colorGroupNames){
			Color color = dataModel.getColorFor(colorGroupValues.get(colorGroupNames.indexOf(groupStr)));

			BallIcon icon = new BallIcon(color);
			JLabel label =  new JLabel(groupStr, icon, SwingConstants.LEFT );
			label.setOpaque(false);
			label.setForeground(super.getForeground());		
			this.add(label, c );

			c.gridy++;			
		}
	}
	
	@Override
	public void setBackground(Color c) {		
		super.setBackground(c);
		for (Component component : this.getComponents()) {
			component.setBackground(c);
		}	
	}
	
	@Override
	public void setForeground(Color c) {
		super.setForeground(c);
		for (Component component : this.getComponents()) {
			component.setForeground(c);
		}
	}
}
