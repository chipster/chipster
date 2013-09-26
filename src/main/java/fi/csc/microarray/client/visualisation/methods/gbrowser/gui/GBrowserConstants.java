package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;

/**
 * Genome browser wide constants.
 *
 */
public class GBrowserConstants {
	
	public static final Color[] charColors = new Color[] {
		new Color(159, 223, 159), // A
		new Color(159, 159, 223), // C
		new Color(191, 191, 191), // G
		new Color(223, 159, 159) // T		
	};
	
	public static final Color BED_COLOR = new Color(1, 119, 183);
	public static final Color COLOR_RED = new Color(0xd53833);
	public static final Color COLOR_BLUE = new Color(0x0177b7);
	public static final Color COLOR_BLUE_BRIGHTER = new Color(0x0199EB);
	public static final Color COLOR_ORANGE = new Color(0xe7881c);
	
	public static final Color SCATTERPLOT_TITLE_COLOR = new Color(0f, 0f, 0f, 0.1f);
	
	public static final Color FORWARD_COLOR = new Color(0,0,0, 64);
	public static final Color REVERSE_COLOR = new Color(
    		GBrowserConstants.COLOR_BLUE.getRed(), 
    		GBrowserConstants.COLOR_BLUE.getGreen(), 
    		GBrowserConstants.COLOR_BLUE.getBlue(),     	    
    		64);

	private static Color coverageColor;
	
	public static Color getCoverageColor() {
		if (coverageColor == null) {
			
			//Calculate the sum of two transparent colors
			BufferedImage img = new BufferedImage(1, 1, BufferedImage.TYPE_4BYTE_ABGR);   		    
		    Graphics imgG = img.getGraphics();
		    
		    imgG.setColor(FORWARD_COLOR);
		    imgG.fillRect (0,0,1,1);		    
		    imgG.setColor(REVERSE_COLOR);
		    imgG.fillRect (0,0,1,1);
		    
		    coverageColor = new Color(img.getRGB(0, 0), true);
		}
		
		return coverageColor;
	}

	// Visibility level thresholds
	public static final int CHANGE_TRACKS_ZOOM_THRESHOLD2 = 10000000;
	public static int SWITCH_VIEWS_AT = 50000;
	public static int SHOW_AVERAGES = 2000;
	public static int SHOW_REFERENCE_AT = 800;
	public static final int SHOW_SNP_AT = 800;
	
	// Read drawing
	public static final int SPACE_BETWEEN_READS = 2;
	public static final int READ_HEIGHT = 4;
	
    // Genome Browser legend
    public static final String GB_LEGEND_A_ICON = 
    		"/gbrowserLegend/a.png";
    public static final String GB_LEGEND_C_ICON = 
    		"/gbrowserLegend/c.png";
    public static final String GB_LEGEND_G_ICON = 
            "/gbrowserLegend/g.png";
    public static final String GB_LEGEND_T_ICON = 
            "/gbrowserLegend/t.png";
    public static final String GB_LEGEND_TOTAL_ICON = 
            "/gbrowserLegend/total.png";
    public static final String GB_LEGEND_FORWARD_ICON = 
            "/gbrowserLegend/forward.png";
    public static final String GB_LEGEND_REVERSE_ICON = 
            "/gbrowserLegend/reverse.png";
    public static final String GB_LEGEND_UTR_ICON = 
            "/gbrowserLegend/utr.png";
    public static final String GB_LEGEND_CDS_ICON = 
            "/gbrowserLegend/cds.png";
    public static final String GB_LEGEND_INTRON_ICON = 
            "/gbrowserLegend/intron.png";
    public static final String GB_LEGEND_INSERTION_ICON = 
            "/gbrowserLegend/insertion.png";
    public static final String GB_LEGEND_DELETION_ICON = 
            "/gbrowserLegend/deletion.png";
    public static final String GB_LEGEND_N_ICON = 
            "/gbrowserLegend/n.png";
    public static final String GB_LEGEND_READ_ICON = 
            "/gbrowserLegend/read.png";
    public static final String GB_LEGEND_MORE_READS_ICON = 
            "/gbrowserLegend/more-reads.png";
    public static final String GB_LEGEND_END_ICON = 
            "/gbrowserLegend/read-end.png";
    public static final String GB_LEGEND_REPEAT_ICON = 
            "/gbrowserLegend/repeat.png";
    public static final String GB_LEGEND_GENE_ICON = 
    		"/gbrowserLegend/gene.png";
}
