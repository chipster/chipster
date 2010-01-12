package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

/**
 * JPanel with black background and rainbow gradient. Values of the different colors
 * are taken from the dataModel and their colors from the color model.
 * 
 * @author Petri KlemelÃ¤
 */
public class ColorGroupsPanel extends JPanel{

	DataModel dataModel;
	Float[] scaleValues;

	public ColorGroupsPanel(DataModel dataModel, List<String> colorGroupList){
		this.dataModel = dataModel;
		this.scaleValues = dataModel.getColorValues();
		//this.setSize(150, 300);
		//this.setPreferredSize(new Dimension(150, 300));
		this.setBackground(Color.black);
		//this.setOpaque(true);				
		this.setLayout(new GridBagLayout());

		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;

		for(String groupStr : colorGroupList){
			Image image = new BufferedImage(16, 16, BufferedImage.TYPE_INT_ARGB);

			Graphics2D g = (Graphics2D)image.getGraphics();

			Color color = dataModel.getColorModel().getColorFor(
					dataModel.convertToScaled(
							dataModel.getColorScaleValues(), 
							colorGroupList.indexOf(groupStr)));

			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON);


			DataPoint.paintBall(2, 2, 12, 12, color, g);

			//				g.setColor(dataModel.getColorModel().getColorFor(
			//						dataModel.convertToScaled(
			//								dataModel.getColorScaleValues(), 
			//								groupValue)));
			//
			//				g.fillRect(4, 4, 10, 10);
			//				
			ImageIcon icon = new ImageIcon(image);			
			JLabel label =  new JLabel(groupStr, icon, SwingConstants.LEFT );
			label.setOpaque(false);
			label.setForeground(Color.white);		
			this.add(label, c );

			c.gridy++;			
		}
	}
}
