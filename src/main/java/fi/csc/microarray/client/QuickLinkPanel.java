package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class QuickLinkPanel extends JPanel implements ActionListener {

	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink emptyLink;
	private JXHyperlink exampleLink;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;
	private JXHyperlink importArrayExpressLink;
	private JXHyperlink importGEOLink;
	private JXHyperlink importUniProtLink;
	private JXHyperlink importEMBLLink;
	private JXHyperlink importTextLink;

	private static final String LINK_WORD = "***";
	
	public QuickLinkPanel() {
		super(new GridBagLayout());

		application = (SwingClientApplication) Session.getSession().getApplication();

		this.setBackground(Color.white);
		
		// Prepare all available links
		exampleLink = createLink("Open example session ");
		importLink = createLink("Import files ");
		importFolderLink = createLink("Import folder ");
		importURLLink = createLink("Import from URL ");
		importArrayExpressLink = createLink("Import from ArrayExpress ");
		importGEOLink = createLink("Import from GEO ");
		importUniProtLink = createLink("Import from UniProt ");
		importEMBLLink = createLink("Import from EMBL ");
		importTextLink = createLink("Create dataset from text ");
		sessionLink = createLink("Open session ");			

		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;

		c.insets.set(5, 10, 5, 10);
		c.gridwidth = 2;	
		this.add(new JLabel("To start working with Chipster, you need to load in data first:"), c);		
		c.gridwidth = 1;
		c.gridy++;

		c.insets.set(0, 10, 0, 0);

		//addLink("To start working with Chipster, you need to load in data first:", (JXHyperlink)null, null, c);
		
		addLink("*** to get familiar with Chipster.", exampleLink, VisualConstants.EXAMPLE_SESSION_ICON, c);

		addLink("*** to continue working on previous sessions.", sessionLink, VisualConstants.OPEN_SESSION_LINK_ICON, c);

		// Common links
		List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
		importLinks.add(importLink);
		importLinks.add(importFolderLink);
		importLinks.add(importURLLink);
		
		// Microarray links
		String linkTemplate = "\n      *** \n      *** \n      *** \n      *** \n      ***";
		if (application.getRequestedModule().equals(ClientApplication.MODULE_MICROARRAY)) {
		    importLinks.add(importArrayExpressLink);
		    importLinks.add(importGEOLink);
		}
		
		// Sequence links
		linkTemplate = "\n      *** \n      *** \n      *** \n      *** \n      *** \n      ***";
		if (application.getRequestedModule().equals(ClientApplication.MODULE_SEQUENCE)) {
	        importLinks.add(importUniProtLink);
	        importLinks.add(importEMBLLink);
		    importLinks.add(importTextLink);
		}

		addLink("Import new data to Chipster: " + linkTemplate, importLinks, VisualConstants.IMPORT_LINK_ICON, c);

		// Panels to take rest of space
		JPanel bottomPanel = new JPanel();
		JPanel rightPanel = new JPanel();
		bottomPanel.setBackground(Color.white);
		rightPanel.setBackground(Color.white);

		c.weightx = 0.0;
		c.weighty = 1.0;
		c.fill = GridBagConstraints.VERTICAL;
		c.gridx = 0;
		c.gridy++;
		this.add(bottomPanel, c);
		c.weightx = 1.0;
		c.weighty = 0.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 3;
		c.gridy = 0;
		this.add(rightPanel, c);

		this.setMinimumSize(new Dimension(0, 0));
		this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));
	}

	private void addLink(String description, JXHyperlink link, ImageIcon icon, GridBagConstraints c) {
		List<JXHyperlink> list = new LinkedList<JXHyperlink>();
		list.add(link);
		addLink(description, list, icon, c);
	}

	private void addLink(String description, List<JXHyperlink> links, ImageIcon icon, GridBagConstraints c) {

		String[] words = description.split(" ");
		int rowChars = 0;
		final int MAX_ROW_CHARS = 40;
		Iterator<JXHyperlink> linkIterator = links.iterator();
		int rowCount = 0;

		c.gridx = 1;
		c.insets.top = 10;
		JPanel row = null;

		for (int i = 0; i < words.length; i++) {
			if (row == null || rowChars + words[i].length() > MAX_ROW_CHARS || words[i].equals("\n")) {

				FlowLayout flow = new FlowLayout(FlowLayout.LEADING);
				flow.setVgap(0);
				flow.setHgap(0);
				row = new JPanel(flow);
				row.setBackground(Color.white);

				c.gridy++;
				this.add(row, c);
				c.insets.top = 0; // After first row

				rowChars = 0;
				rowCount++;
			}

			if (words[i].equals(LINK_WORD)) {

				JXHyperlink link = linkIterator.next();
				rowChars += link.getText().length() + 1;
				row.add(link);
				// row.add(new JLabel(link.getText()));
			} else if (!words[i].equals("\n")) {
				JLabel text = new JLabel(words[i] + " ");
				row.add(text);
				rowChars += words[i].length() + 1;
			}
		}

		c.gridy -= (rowCount - 1);
		c.gridheight = rowCount;
		c.gridx = 0;
		c.insets.top = 15;
		this.add(new JLabel(icon), c);
		c.gridy += (rowCount - 1);
		c.gridheight = 1;
	}

	private JXHyperlink createLink(String text){
		JXHyperlink link = new JXHyperlink();
		link.setText(text);
		link.addActionListener(this);
		link.setBorder(null);
		link.setMargin(new Insets(0, 0, 0, 0));
		return link;
	}

	public void actionPerformed(ActionEvent e) {
		try {
			if (e.getSource() == sessionLink) {
				application.loadSession();
			} else if (e.getSource() == importLink) {
				application.openFileImport();
			} else if (e.getSource() == importURLLink) {
				application.openURLImport();
			} else if (e.getSource() == importFolderLink) {
				application.openDirectoryImportDialog();
			} else if (e.getSource() == importArrayExpressLink) {
				Operation importOperation = new Operation(application.locateOperationDefinition(MicroarrayModule.IMPORT_CAT, MicroarrayModule.IMPORT_FROM_ARRAYEXPRESS_NAME), new DataBean[] {});
				application.openDatabaseImport("ArrayExpress", importOperation);
			} else if (e.getSource() == importGEOLink) {
				Operation importOperation = new Operation(application.locateOperationDefinition(MicroarrayModule.IMPORT_CAT, MicroarrayModule.IMPORT_FROM_GEO_NAME), new DataBean[] {});
				application.openDatabaseImport("GEO", importOperation);
			} else if (e.getSource() == importTextLink) {
			    application.openCreateFromTextDialog();
			} else if (e.getSource() == emptyLink) {

			} else if (e.getSource() == exampleLink) {
				URL url = new URL("https://extras.csc.fi/biosciences/Chipster/sessionIlluminaHuman6v1teratospermia.cs");
				application.loadSessionFrom(url);
			}
		} catch (Exception ex) {
			application.reportException(ex);
		}
	}
}
