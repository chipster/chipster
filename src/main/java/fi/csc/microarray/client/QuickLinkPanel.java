package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.util.LinkedList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.util.LinkUtil;
import fi.csc.microarray.util.Strings;

@SuppressWarnings("serial")
public class QuickLinkPanel extends JPanel {
	
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(QuickLinkPanel.class); 

	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink localSessionLink;
	private JXHyperlink importLink;
	private JXHyperlink exampleLink;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;
	private JXHyperlink importURLToServerLink;

	public QuickLinkPanel() {
		super(new GridBagLayout());

		application = (SwingClientApplication) Session.getSession().getApplication();

		this.setBackground(Color.white);

		//
		// Prepare all available links
		//
		
		exampleLink = LinkUtil.createLink("Open example session", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.loadSession(true, true, true);							
			}
		});
		
		importLink = LinkUtil.createLink("Import files", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openFileImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		importFolderLink = LinkUtil.createLink("Import folder", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.openDirectoryImportDialog();
			}
		});
		importURLLink = LinkUtil.createLink("Import from URL to client", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openURLImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		
		importURLToServerLink = LinkUtil.createLink("Import from URL directly to server", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				BasicModule.importFromUrlToServer();
			}
		});
		
		sessionLink = LinkUtil.createLink("open cloud session", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.loadSession(true);
			}
		});
		localSessionLink = LinkUtil.createLink("Open local session", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.loadSession(false);
			}
		});

		
		
		// Draw panel
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;

		c.insets.set(5, 10, 5, 10);
		c.gridwidth = 2;
		this.add(new JLabel("To start working with " + Session.getSession().getPrimaryModule().getDisplayName() + ", you need to load in data first."), c);
		c.gridwidth = 1;
		c.gridy++;

		c.insets.set(0, 10, 0, 0);
		
		addLink("*** to get familiar with " + Session.getSession().getPrimaryModule().getDisplayName(), exampleLink, VisualConstants.EXAMPLE_SESSION_ICON, c, this);
	
		String cloudSessionsString = "*** to continue working on previous sessions. You can also *** from the server.";
		String localSessionsString = "*** to continue working on previous sessions.";
		
		if (application.getSessionManager().areCloudSessionsEnabled()) {
			List<JXHyperlink> openLinks = new LinkedList<JXHyperlink>();
			openLinks.add(localSessionLink);
			openLinks.add(sessionLink);		
			addLinks(cloudSessionsString, openLinks, VisualConstants.OPEN_SESSION_LINK_ICON, c, this);
		} else {
			List<JXHyperlink> openLinks = new LinkedList<JXHyperlink>();
			openLinks.add(localSessionLink);
			addLinks(localSessionsString, openLinks, VisualConstants.OPEN_SESSION_LINK_ICON, c, this);	
		}

		
		// common links
		List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
		importLinks.add(importLink);
		importLinks.add(importFolderLink);
		importLinks.add(importURLLink);
		importLinks.add(importURLToServerLink);

		// module specific links
		if (!application.isStandalone()) {
			Module primaryModule = Session.getSession().getPrimaryModule();
			primaryModule.addImportLinks(this, importLinks);
		}
		
		String linkTemplate = Strings.repeat("\n      *** ", importLinks.size());
		addLinks("Import new data to " + Session.getSession().getPrimaryModule().getDisplayName() + ": " + linkTemplate, importLinks, VisualConstants.IMPORT_LINK_ICON, c, this);

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
	
	
	private final static int MAX_ROW_CHARS = 53;

	public static void addLink(String description, JXHyperlink link, String iconPath, GridBagConstraints c, JComponent component) {
	
		LinkUtil.addLink(description, link, VisualConstants.getIcon(iconPath), c, component, MAX_ROW_CHARS, Color.white);
	}

	
	public static void addLinks(String description, List<JXHyperlink> links, String iconPath, GridBagConstraints c, JComponent component) {
		LinkUtil.addLinks(description, links, VisualConstants.getIcon(iconPath), c, component, MAX_ROW_CHARS, Color.white);
	}

}
