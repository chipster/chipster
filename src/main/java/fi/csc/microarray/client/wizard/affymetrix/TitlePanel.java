package fi.csc.microarray.client.wizard.affymetrix;

import java.awt.*;
import javax.swing.UIManager;
import javax.swing.*;
import javax.swing.border.*;


public class TitlePanel extends JPanel {
	private JSeparator separator;

	private JLabel textLabel;

	private JPanel titlePanel;

	private JPanel secondaryPanel;

	public TitlePanel() {
		super();

		//Color grayColor = Color.gray;
		Color backGroundColor = UIManager.getColor("SimpleInternalFrame.activeTitleBackground");
		Color titleFontColor = UIManager.getColor("SimpleInternalFrame.activeTitleForeground");
		
		titlePanel = new javax.swing.JPanel();
		textLabel = new javax.swing.JLabel();
		separator = new javax.swing.JSeparator();

		setLayout(new java.awt.BorderLayout());

		titlePanel.setLayout(new java.awt.BorderLayout());
		titlePanel.setBackground(backGroundColor);
		//titlePanel.setBackground(grayColor.brighter());

		//textLabel.setBackground(grayColor.brighter());
		textLabel.setBackground(backGroundColor);
		textLabel.setForeground(titleFontColor);
		textLabel.setFont(new Font("MS Sans Serif", Font.BOLD, 14));
		textLabel.setBorder(new EmptyBorder(new Insets(3, 5, 3, 3)));
		textLabel.setOpaque(true);

		titlePanel.add(textLabel, BorderLayout.CENTER);
		titlePanel.add(separator, BorderLayout.SOUTH);

		add(titlePanel, BorderLayout.NORTH);

		secondaryPanel = new JPanel();
		add(secondaryPanel, BorderLayout.WEST);
	}

	public void setTitle(String title) {
		textLabel.setText(title);
	}

	public void addContent(JPanel contentPanel) {
		secondaryPanel.add(contentPanel, BorderLayout.NORTH);
	}
}
