package fi.csc.microarray.util;

import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fi.csc.microarray.client.dataimport.ImportUtils;

public class ImageExportUtils {
	
	enum Resolution {
		RES1("1X", 1),
		RES2("2X", 2),
		RES4("4X", 4),
		RES8("8X", 8);			
		
		private String name;
		private int resolution;
		
		Resolution(String name, int resolution) {
			this.name = name;
			this.resolution = resolution;
		}
		
		public String toString() {
			return name;
		}
	}
	
	public static class ResolutionAccessory extends JPanel {
		
		private JComboBox<Resolution> comboBox;

		public ResolutionAccessory() {
			JLabel label = new JLabel("Resolution: ");
			comboBox = new JComboBox<>(Resolution.values());			

			this.add(label);
			this.add(comboBox);
		}
		
		public int getResolution() {
			Resolution res = (Resolution) comboBox.getSelectedItem();
			return res.resolution;
		}
	}
	public static BufferedImage componentToImage(Component component, int resolution) throws IOException
	{
	    BufferedImage img = new BufferedImage(component.getWidth() * resolution, component.getHeight() * resolution, BufferedImage.TYPE_INT_ARGB_PRE);
	    Graphics2D g2 = (Graphics2D) img.getGraphics();
	    g2.scale(resolution, resolution);
	    g2.setColor(component.getForeground());
	    g2.setFont(component.getFont());
	    component.paintAll(g2);
	    return img;
	}
	
	public static JFileChooser getSaveFileChooser() {
		
		JFileChooser fileChooser = ImportUtils.getFixedFileChooser();
		
		String[] extensions = { "png" };
		fileChooser.setFileFilter(new GeneralFileFilter("PNG Image Files", extensions));
		fileChooser.setSelectedFile(new File("genome-browser.png"));
		fileChooser.setAcceptAllFileFilterUsed(false);
		fileChooser.setMultiSelectionEnabled(false);		
		fileChooser.setApproveButtonText("Save");
		fileChooser.setAccessory(new ResolutionAccessory());
		
		return fileChooser;
	}

	public static JFileChooser saveComponent(Component component, JFileChooser fileChooser) throws IOException {
		
		if (fileChooser == null) {
			fileChooser = ImageExportUtils.getSaveFileChooser();
		}
		
		int option = fileChooser.showOpenDialog(component);
		if (option == JFileChooser.APPROVE_OPTION) {					
			ResolutionAccessory accessory =  (ResolutionAccessory) fileChooser.getAccessory();
			
			BufferedImage img = ImageExportUtils.componentToImage(component, accessory.getResolution());
			
			ImageIO.write(img, "png", fileChooser.getSelectedFile());
		}
				
		
		return fileChooser;
	}

	public static int printComponent(Graphics g, PageFormat pf, int page, Component component) {
	    // We have only one page, and 'page'
	    // is zero-based
	    if (page > 0) {
	         return Printable.NO_SUCH_PAGE;
	    }
	    
	    Graphics2D g2d = (Graphics2D)g;


	    // User (0,0) is typically outside the
	    // imageable area, so we must translate
	    // by the X and Y values in the PageFormat
	    // to avoid clipping.
	    g2d.translate(pf.getImageableX(), pf.getImageableY());
	    
	    //Scale component image to one page size
	    double scaleX = pf.getImageableWidth() / component.getWidth();
	    double scaleY = pf.getImageableHeight() / component.getHeight();
	    double scale = Math.min(scaleX, scaleY);
	    
	    g2d.scale(scale, scale);

	    component.paintAll(g);

	    // tell the caller that this page is part
	    // of the printed document
	    return Printable.PAGE_EXISTS;
	}

	public static void printComponent(Printable component) throws PrinterException {
		PrinterJob job = PrinterJob.getPrinterJob();
		job.setPrintable(component);
		
		boolean doPrint = job.printDialog();
		if (doPrint) {
			job.print();
		}
	}
}
