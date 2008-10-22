package fi.csc.microarray.wizard.affymetrix;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;


public class NormalisationPanel extends TitlePanel {
 
	private JComboBox methodComboBox = new JComboBox();
	private SelectionSet selectionSet = new SelectionSet();
	private String methodName = null;

	public NormalisationPanel() {
		super();

		setTitle("Wizard: Normalisation");

		JPanel contentPanel = getContentPanel();
		addContent(contentPanel);
	}

	private JPanel getContentPanel() {
		selectionSet.add("MAS5", "mas5");
		selectionSet.add("RMA", "rma");
		selectionSet.add("GCRMA", "gcrma");
		selectionSet.initComboBox(methodComboBox);
		
		JPanel cPanel = new JPanel();
		cPanel.setLayout(new BorderLayout());

		JPanel upperPanel = new JPanel();
		upperPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		upperPanel.setBorder(new EmptyBorder(new Insets(50, 25, 10, 0)));
		cPanel.add(upperPanel, BorderLayout.CENTER);

		upperPanel.add(new JLabel("Normalisation method:"));
		upperPanel.add(methodComboBox);
		methodComboBox.setEditable(false);
		methodComboBox.addActionListener(new MethodSelector());
		methodComboBox.setSelectedIndex(1);
		methodName = "rma"; // default value

		JPanel lowerPanel = new JPanel();
		cPanel.add(lowerPanel, BorderLayout.SOUTH);
		lowerPanel.setLayout(new BoxLayout(lowerPanel, BoxLayout.Y_AXIS));
		lowerPanel.setBorder(new EmptyBorder(new Insets(10, 25, 10, 10)));
		
		lowerPanel.add(new JLabel("Before any further analysis, your data needs to be normalized. There are"));
		lowerPanel.add(new JLabel("three different methods you can choose from: MAS5, RMA and GCRMA. MAS5 uses"));
		lowerPanel.add(new JLabel("both perfect match and mismatch probes for estimating the gene expression."));
		lowerPanel.add(new JLabel("RMA and GCRMA disregard mismatch probes while estimating the gene"));
		lowerPanel.add(new JLabel("expression. In addition, GCRMA takes the GC-content of the probes into"));
		lowerPanel.add(new JLabel("account. RMA is the typical choice."));
		
		return cPanel;
	}

	public String getSelected() {
		return methodName;
	}

	class MethodSelector implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			Object obj = methodComboBox.getSelectedItem();
			methodName = selectionSet.getInternalValue(obj.toString());
		}
	}
}
