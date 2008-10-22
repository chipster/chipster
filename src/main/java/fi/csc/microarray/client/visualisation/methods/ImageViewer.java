package fi.csc.microarray.client.visualisation.methods;

import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;

public class ImageViewer extends Visualisation implements MouseListener {
	private static final Logger logger = Logger.getLogger(ImageViewer.class);
	
	private class ImagePanel extends JPanel{
		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			Dimension size = ImageViewer.this.calculateSize();
			
			//Has to be done after the visualisation has its parent to know the size
			ImageViewer.this.updateCursor();
			
			Graphics2D g2 = (Graphics2D)g;
			g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
					RenderingHints.VALUE_INTERPOLATION_BILINEAR);
			g2.drawImage(original, 0, 0, (int)size.getWidth(), (int)size.getHeight(), null);
		}
	}
	
	private Image original;
	private JScrollPane scroller;
	private boolean isScaledMode = true;
	private ImagePanel imagePanel = new ImagePanel();
	
	public ImageViewer(VisualisationFrame frame){
		super(frame);
		imagePanel.addMouseListener(this);
	}
	
	private static final Cursor ZOOM_IN_CURSOR = Toolkit.getDefaultToolkit().
		createCustomCursor(VisualConstants.ZOOM_IN_CURSOR_IMAGE.getImage(), new Point(5,2), "ZoomIn");
	private static final Cursor ZOOM_OUT_CURSOR = Toolkit.getDefaultToolkit().
		createCustomCursor(VisualConstants.ZOOM_OUT_CURSOR_IMAGE.getImage(), new Point(5,2), "ZoomOut");

	private double calculateScale(){
						
		Dimension targetSize = scroller.getViewport().getExtentSize();
		double xScale = targetSize.getWidth() / original.getWidth(null);
		double yScale = targetSize.getHeight() / original.getHeight(null);
		
	//	logger.debug("Original x: " + original.getWidth(null));
	//	logger.debug("Target x  : " + targetSize.getWidth());
		
		return xScale < yScale ?  xScale : yScale; 
	}
	
	private Dimension calculateSize(){
		double scale = this.calculateScale();
		if(scale >= 1.0 || !isScaledMode){
			scale = 1.0;
			isScaledMode = false;
		}
		
		return new Dimension((int)(original.getWidth(null)*scale), 
				(int)(original.getHeight(null)*scale));
	}
	
	private void updateCursor(){
		if(calculateScale() < 1.0){			
			if(isScaledMode){
				logger.debug("Zoom in cursor");
				imagePanel.setCursor(ZOOM_IN_CURSOR);
			} else {
				logger.debug("Zoom out cursor");
				imagePanel.setCursor(ZOOM_OUT_CURSOR);
			}
		} else {
			isScaledMode = false;
			logger.debug("Null cursor");
			imagePanel.setCursor(null);
		}
	}

	public void mouseClicked(MouseEvent e) {
		logger.debug("mouse clicked: " + e);
		if(calculateScale() < 1.0){
			isScaledMode = !isScaledMode;		
		}
		this.updateCursor();
				
		//Both are needed
		imagePanel.setPreferredSize(this.calculateSize());
		imagePanel.setSize(this.calculateSize());
		
		//logger.debug("imagePanel preferredSize: " + imagePanel.getPreferredSize());
		//logger.debug("imagePanel size: " + imagePanel.getSize());
		logger.debug("isScaledMode: " + isScaledMode);
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
	}



	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		byte[] bytes = data.getContents();
		if (bytes != null) {
			this.original = new ImageIcon(bytes).getImage();
			this.scroller = new JScrollPane();			
			scroller.setViewportView(imagePanel);
			scroller.setBorder(BorderFactory.createEmptyBorder());
			isScaledMode = true;
			return scroller;
		}
		return this.getDefaultVisualisation();
	}
	
	

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("image/jpeg", "image/png", "image/gif");
	}
}