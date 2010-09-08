package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.util.Strings;

@SuppressWarnings("serial")
public class QuickLinkPanel extends JPanel {

	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink exampleLink;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;

	private static final String LINK_WORD = "***";

	public QuickLinkPanel() {
		super(new GridBagLayout());

		application = (SwingClientApplication) Session.getSession().getApplication();

		this.setBackground(Color.white);

		//
		// Prepare all available links
		//
		
		// Check if example session is available
		exampleLink = null;
		try {
			final URL url = Session.getSession().getModules().getPrimaryModule().getExampleSessionUrl();
			if (url != null) {
				exampleLink = createLink("Open example session ", new AbstractAction() {
					@Override
					public void actionPerformed(ActionEvent e) {
						try {
							application.loadSessionFrom(url);
							hackExampleSessionWithTypeTags();
						} catch (Exception exception) {
							application.reportException(exception);
						}
					}
				});
			}
		} catch (MalformedURLException mue) {
			// ignore and let exampleLink be null
		}

		importLink = createLink("Import files ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openFileImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		importFolderLink = createLink("Import folder ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.openDirectoryImportDialog();
			}
		});
		importURLLink = createLink("Import from URL ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openURLImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		sessionLink = createLink("Open session ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.loadSession();
			}
		});

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

		if (exampleLink != null) {
			addLink("*** to get familiar with Chipster.", exampleLink, VisualConstants.EXAMPLE_SESSION_ICON, c);
		}
		
		addLink("*** to continue working on previous sessions.", sessionLink, VisualConstants.OPEN_SESSION_LINK_ICON, c);

		// common links
		List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
		importLinks.add(importLink);
		importLinks.add(importFolderLink);
		importLinks.add(importURLLink);

		// module specific links
		Module primaryModule = Session.getSession().getModules().getPrimaryModule();
		primaryModule.addImportLinks(this, importLinks);
		
		String linkTemplate = Strings.repeat("\n      *** ", importLinks.size());
		addLinks("Import new data to Chipster: " + linkTemplate, importLinks, VisualConstants.IMPORT_LINK_ICON, c);

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
		addLinks(description, list, icon, c);
	}

	private void addLinks(String description, List<JXHyperlink> links, ImageIcon icon, GridBagConstraints c) {

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

	public JXHyperlink createLink(String text, Action action) {
		JXHyperlink link = new JXHyperlink();
		link.setBorder(null);
		link.setMargin(new Insets(0, 0, 0, 0));
		link.setAction(action);
		link.setText(text); // must be after setAction
		return link;
	}
	
	private void hackExampleSessionWithTypeTags() {
		DataFolder folder = application.getDataManager().getRootFolder().getChildFolder("IlluminaTeratospermiaHuman6v1_BS1");
		
		for (DataItem item : folder.getChildren()) {
			if (item instanceof DataBean) {
				((DataBean)item).addTypeTag(BasicModule.TABLE_WITHOUT_HEADER_TT);
			}
		}
		
	}

}
