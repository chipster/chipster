package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;

public class GBrowserLegend extends JPanel {
	
	private GBrowser browser;

	public GBrowserLegend(GBrowser browser) {
		
		this.browser = browser;
		
		this.setLayout(new MigLayout());
		this.setBackground(Color.white);
		
		this.add(GBrowserSettings.createTitle("Annotation track"), "span 2, wrap");	
		
		addIcon(GBrowserConstants.GB_LEGEND_GENE_ICON, "Gene");
		addIcon(GBrowserConstants.GB_LEGEND_CDS_ICON, "Coding sequence");
		addIcon(GBrowserConstants.GB_LEGEND_UTR_ICON, "Untranslated region");
		addIcon(GBrowserConstants.GB_LEGEND_INTRON_ICON, "Intron");
		addIcon(GBrowserConstants.GB_LEGEND_END_ICON, "Transcript end");
		addIcon(GBrowserConstants.GB_LEGEND_REPEAT_ICON, "Low complexity region");
		
		this.add(GBrowserSettings.createTitle("Sample tracks"), "span 2, wrap");
		
		addIcon(GBrowserConstants.GB_LEGEND_READ_ICON, "Read");
		addIcon(GBrowserConstants.GB_LEGEND_MULTIMAPPING_ICON, "Multimapping read");
		addIcon(GBrowserConstants.GB_LEGEND_MORE_READS_ICON, "More reads to show");
		addIcon(GBrowserConstants.GB_LEGEND_INSERTION_ICON, "Insertion");
		addIcon(GBrowserConstants.GB_LEGEND_DELETION_ICON, "Deletion");
		
		addIcon(GBrowserConstants.GB_LEGEND_A_ICON, "Adenine");
		addIcon(GBrowserConstants.GB_LEGEND_C_ICON, "Cytosine");
		addIcon(GBrowserConstants.GB_LEGEND_G_ICON, "Guanine");
		addIcon(GBrowserConstants.GB_LEGEND_T_ICON, "Thymine");
		addIcon(GBrowserConstants.GB_LEGEND_N_ICON, "N");
		
		addIcon(GBrowserConstants.GB_LEGEND_TOTAL_ICON, "Total coverage");
		addIcon(GBrowserConstants.GB_LEGEND_FORWARD_ICON, "Forward coverage");
		addIcon(GBrowserConstants.GB_LEGEND_REVERSE_ICON, "Reverse coverage");	
	}

	private void addIcon(String path, String text) {

		ImageIcon icon = browser.getIcon(path);
		
		this.add(new JLabel(icon));
		this.add(new JLabel(text), "wrap");
	}
}
