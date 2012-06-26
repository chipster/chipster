package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;

import javax.swing.ImageIcon;

import fi.csc.microarray.constants.VisualConstants;

/**
 * Genome browser wide constants.
 *
 */
public class GenomeBrowserConstants {
	
	public static final Color[] charColors = new Color[] {
		new Color(159, 223, 159), // A
		new Color(159, 159, 223), // C
		new Color(191, 191, 191), // G
		new Color(223, 159, 159) // T
		
//		new Color(64, 192, 64, 128), // A
//		new Color(64, 64, 192, 128), // C
//		new Color(128, 128, 128, 128), // G
//		new Color(192, 64, 64, 128) // T
	};

	// Visibility level thresholds
	public static final int CHANGE_TRACKS_ZOOM_THRESHOLD2 = 10000000;
	public static int SWITCH_VIEWS_AT = 50000;
	public static int SHOW_REFERENCE_AT = 800;
	public static final int SHOW_SNP_AT = 800;
	
	// Read drawing
	public static final int SPACE_BETWEEN_READS = 2;
	public static final int READ_HEIGHT = 4;
	
    // Genome Browser legend
    public static final ImageIcon GB_LEGEND_A_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/a.png"));
    public static final ImageIcon GB_LEGEND_C_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/c.png"));
    public static final ImageIcon GB_LEGEND_G_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/g.png"));
    public static final ImageIcon GB_LEGEND_T_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/t.png"));
    public static final ImageIcon GB_LEGEND_TOTAL_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/total.png"));
    public static final ImageIcon GB_LEGEND_FORWARD_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/forward.png"));
    public static final ImageIcon GB_LEGEND_REVERSE_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/reverse.png"));
    public static final ImageIcon GB_LEGEND_UTR_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/utr.png"));
    public static final ImageIcon GB_LEGEND_CDS_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/cds.png"));
    public static final ImageIcon GB_LEGEND_INTRON_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/intron.png"));
    public static final ImageIcon GB_LEGEND_INSERTION_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/insertion.png"));
    public static final ImageIcon GB_LEGEND_DELETION_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/deletion.png"));
    public static final ImageIcon GB_LEGEND_N_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/n.png"));
    public static final ImageIcon GB_LEGEND_READ_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/read.png"));
    public static final ImageIcon GB_LEGEND_END_ICON = 
            new ImageIcon(VisualConstants.class.getResource("/gbrowserLegend/read-end.png"));


}
