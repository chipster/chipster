package fi.csc.microarray.wizard.affymetrix;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;


public class GroupTestPanel extends TitlePanel {
 
	private JComboBox testComboBox = new JComboBox();
	private SelectionSet selectionSet = new SelectionSet();
	private String testName = null;

	public GroupTestPanel() {
		super();

		setTitle("Wizard: Group test");

		JPanel contentPanel = getContentPanel();
		addContent(contentPanel);
	}

	private JPanel getContentPanel() {
		
		selectionSet.add("Empirical Bayes", "empiricalBayes+BH");
		selectionSet.add("ANOVA", "ANOVA+BH");
		selectionSet.add("Kruskal-Wallis", "Kruskal-Wallis+BH");
		selectionSet.initComboBox(testComboBox);
		
		JPanel cPanel = new JPanel();
		cPanel.setLayout(new BorderLayout());
		
		JPanel upperPanel = new JPanel();
		upperPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		upperPanel.setBorder(new EmptyBorder(new Insets(50, 25, 10, 0)));
		cPanel.add(upperPanel, BorderLayout.CENTER);
		
		upperPanel.add(new JLabel("Statistical test:"));
		upperPanel.add(testComboBox);
		testComboBox.setEditable(false);
		testComboBox.addActionListener(new TestSelector());

		testName = "empiricalBayes+BH"; // default value

		JPanel lowerPanel = new JPanel();
		cPanel.add(lowerPanel, BorderLayout.SOUTH);
		lowerPanel.setLayout(new BoxLayout(lowerPanel, BoxLayout.Y_AXIS));
		lowerPanel.setBorder(new EmptyBorder(new Insets(10, 25, 10, 10)));
		
		lowerPanel.add(new JLabel("To find the differentially expressed genes, a statistical test needs to be"));
		lowerPanel.add(new JLabel("applied. There are three different approaches to testing: empirical Bayes,"));
		lowerPanel.add(new JLabel("analysis of variance (ANOVA), and Kruskall-Wallis. The first two assume that"));
		lowerPanel.add(new JLabel("gene expression values are approximately normally distributed. The last does"));
		lowerPanel.add(new JLabel("not assume normality. Empirical Bayes is the typical choice at this stage of"));
		lowerPanel.add(new JLabel("an analysis."));
		
		return cPanel;
	}

	public String getSelected() {
		return testName;
	}

	class TestSelector implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			Object obj = testComboBox.getSelectedItem();
			testName = selectionSet.getInternalValue(obj.toString());
		}
	}
}
