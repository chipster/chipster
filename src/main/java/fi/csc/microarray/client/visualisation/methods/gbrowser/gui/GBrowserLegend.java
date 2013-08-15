package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;

public class GBrowserLegend extends JPanel {
	
	private JPanel panel;
	private GBrowser browser;

	public GBrowserLegend(GBrowser browser) {
		
		this.browser = browser;
		
		panel = new JPanel();
		panel.setBackground(Color.white);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		//c.gridheight = 32;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets = new Insets(7, 7, 7, 7);
		c.weighty = 0;
		
		addIcon(GBrowserConstants.GB_LEGEND_GENE_ICON, "Gene", c);
		addIcon(GBrowserConstants.GB_LEGEND_CDS_ICON, "Coding sequence", c);
		addIcon(GBrowserConstants.GB_LEGEND_UTR_ICON, "Untranslated region", c);
		addIcon(GBrowserConstants.GB_LEGEND_INTRON_ICON, "Intron", c);
		addIcon(GBrowserConstants.GB_LEGEND_END_ICON, "Transcript end", c);
		addIcon(GBrowserConstants.GB_LEGEND_REPEAT_ICON, "Low complexity region", c);
		
		addIcon(GBrowserConstants.GB_LEGEND_READ_ICON, "Read", c);
		addIcon(GBrowserConstants.GB_LEGEND_MORE_READS_ICON, "More reads to show", c);
		addIcon(GBrowserConstants.GB_LEGEND_INSERTION_ICON, "Insertion", c);
		addIcon(GBrowserConstants.GB_LEGEND_DELETION_ICON, "Deletion", c);
		
		addIcon(GBrowserConstants.GB_LEGEND_A_ICON, "Adenine", c);
		addIcon(GBrowserConstants.GB_LEGEND_C_ICON, "Cytosine", c);
		addIcon(GBrowserConstants.GB_LEGEND_G_ICON, "Guanine", c);
		addIcon(GBrowserConstants.GB_LEGEND_T_ICON, "Thymine", c);
		addIcon(GBrowserConstants.GB_LEGEND_N_ICON, "N", c);
		
		addIcon(GBrowserConstants.GB_LEGEND_TOTAL_ICON, "Total coverage", c);
		addIcon(GBrowserConstants.GB_LEGEND_FORWARD_ICON, "Forward coverage", c);
		addIcon(GBrowserConstants.GB_LEGEND_REVERSE_ICON, "Reverse coverage", c);
		
		c.fill = GridBagConstraints.VERTICAL;
		c.weighty = 1.0;
		
		JPanel spaceFiller = new JPanel();
		spaceFiller.setBackground(Color.white);
		
		panel.add(spaceFiller, c);

		JScrollPane scroller = new JScrollPane(panel);
		scroller.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		this.setLayout(new BorderLayout());
		this.add(scroller, BorderLayout.CENTER);
		
	}

	private void addIcon(String path, String text, GridBagConstraints c) {

		ImageIcon icon = browser.getIcon(path);
		
		c.gridx = 0;
		panel.add(new JLabel(icon), c);
		
		c.gridx = 1;
		panel.add(new JLabel(text), c);
		
		c.gridy++;
	}
}
