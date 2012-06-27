package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fi.csc.microarray.constants.VisualConstants;

public class GBrowserLegend extends JPanel {

	public GBrowserLegend() {
		
		this.setBackground(Color.white);
		this.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		//c.gridheight = 32;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets = new Insets(10, 10, 10, 10);
		c.weighty = 0;
		
		addIcon(VisualConstants.GB_LEGEND_UTR_ICON, "Untranslated region", c);
		addIcon(VisualConstants.GB_LEGEND_START_CODON_ICON, "Start codon", c);
		addIcon(VisualConstants.GB_LEGEND_CDS_ICON, "Coding sequence", c);
		addIcon(VisualConstants.GB_LEGEND_INTRON_ICON, "Intron", c);
		
		addIcon(VisualConstants.GB_LEGEND_READ_ICON, "Read", c);
		addIcon(VisualConstants.GB_LEGEND_END_ICON, "3' end", c);
		addIcon(VisualConstants.GB_LEGEND_INSERTION_ICON, "Insertion", c);
		addIcon(VisualConstants.GB_LEGEND_DELETION_ICON, "Deletion", c);
		
		addIcon(VisualConstants.GB_LEGEND_A_ICON, "Adenosine", c);
		addIcon(VisualConstants.GB_LEGEND_C_ICON, "Cytidine", c);
		addIcon(VisualConstants.GB_LEGEND_G_ICON, "Guanosine", c);
		addIcon(VisualConstants.GB_LEGEND_T_ICON, "Thymidine", c);
		addIcon(VisualConstants.GB_LEGEND_N_ICON, "N", c);
		
		addIcon(VisualConstants.GB_LEGEND_FORWARD_ICON, "Forward coverage", c);
		addIcon(VisualConstants.GB_LEGEND_REVERSE_ICON, "Reverse coverage", c);
		c.fill = GridBagConstraints.VERTICAL;
		c.weighty = 1.0;
		
		JPanel spaceFiller = new JPanel();
		spaceFiller.setBackground(Color.white);
		
		this.add(spaceFiller, c);
		
	}

	private void addIcon(ImageIcon icon, String text, GridBagConstraints c) {

		c.gridx = 0;
		this.add(new JLabel(icon), c);
		
		c.gridx = 1;
		this.add(new JLabel(text), c);
		
		c.gridy++;
	}
}
