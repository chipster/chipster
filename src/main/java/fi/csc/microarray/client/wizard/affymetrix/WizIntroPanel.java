package fi.csc.microarray.client.wizard.affymetrix;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;


public class WizIntroPanel extends TitlePanel {

	public WizIntroPanel() {
		super();

		setTitle("Wizard: Introduction");

		JPanel contentPanel = getContentPanel();
		contentPanel.setBorder(new EmptyBorder(new Insets(50, 15, 10, 10)));

		addContent(contentPanel);
	}

	private JPanel getContentPanel() {
		JPanel cPanel = new JPanel();
		cPanel.setLayout(new BoxLayout(cPanel, BoxLayout.Y_AXIS));

		JPanel textPane = new JPanel();
		textPane.setLayout(new java.awt.GridLayout(0, 1));
		cPanel.add(textPane);

		textPane.add(new JLabel("This wizard enables you to find the differentially expressed genes from"));
		textPane.add(new JLabel("your Affymetrix data with just a few mouse clicks. This wizard will put your"));
		textPane.add(new JLabel("CEL-files through an analysis that automatically normalizes the data, finds"));
		textPane.add(new JLabel("the differentially expressed genes between two or more treatments or groups,"));
		textPane.add(new JLabel("and performs hierarchical clustering for the data."));

		return cPanel;
	}
}
