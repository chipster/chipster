package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.util.Strings;

@SuppressWarnings("serial")
public class QuickLinkPanel extends JPanel {

	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink exampleLink;
	private JXHyperlink exampleLinkAlternative;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;

	public QuickLinkPanel() {
		super(new GridBagLayout());

		application = (SwingClientApplication) Session.getSession().getApplication();

		this.setBackground(Color.white);

		//
		// Prepare all available links
		//
		
		// Check if example session is available
		exampleLink = null;
		exampleLinkAlternative = null;

		try {
			final URL[] urls = Session.getSession().getPrimaryModule().getExampleSessionUrls(application.isStandalone);
			if (urls != null) {

				if (urls.length == 1) {
					exampleLink = LinkUtil.createLink("Example session ", new AbstractAction() {
						@Override
						public void actionPerformed(ActionEvent e) {
							try {
								application.loadSessionFrom(urls[0]);
							} catch (Exception exception) {
								application.reportException(exception);
							}
						}
					});
				}
				
				if (urls.length == 2) {
					exampleLink = LinkUtil.createLink("microarray", new AbstractAction() {
						@Override
						public void actionPerformed(ActionEvent e) {
							try {
								application.loadSessionFrom(urls[0]);
							} catch (Exception exception) {
								application.reportException(exception);
							}
						}
					});
					
					exampleLinkAlternative = LinkUtil.createLink("NGS", new AbstractAction() {
						@Override
						public void actionPerformed(ActionEvent e) {
							try {
								application.loadSessionFrom(urls[1]);
							} catch (Exception exception) {
								application.reportException(exception);
							}
						}
					});
				}
			}
		} catch (MalformedURLException mue) {
			// ignore and let exampleLink be null
		}
		
		importLink = LinkUtil.createLink("Import files ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openFileImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		importFolderLink = LinkUtil.createLink("Import folder ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.openDirectoryImportDialog();
			}
		});
		importURLLink = LinkUtil.createLink("Import from URL ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					application.openURLImport();
				} catch (Exception exception) {
					application.reportException(exception);
				}
			}
		});
		sessionLink = LinkUtil.createLink("Open session ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				application.loadSession();
			}
		});

		
		
		// Draw panel
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;

		c.insets.set(5, 10, 5, 10);
		c.gridwidth = 2;
		this.add(new JLabel("To start working with " + Session.getSession().getPrimaryModule().getDisplayName() + ", you need to load in data first:"), c);
		c.gridwidth = 1;
		c.gridy++;

		c.insets.set(0, 10, 0, 0);

		if (exampleLink != null) {
						
			if (exampleLinkAlternative == null) {
				
				exampleLink.setText("Open example session ");
				addLink("*** to get familiar with " + Session.getSession().getPrimaryModule().getDisplayName() + ". " , exampleLink, VisualConstants.EXAMPLE_SESSION_ICON, c, this);

			} else {

				List<JXHyperlink> exampleLinks = new LinkedList<JXHyperlink>();
				exampleLinks.add(exampleLink);
				exampleLinks.add(exampleLinkAlternative);

				addLinks("Open example session (*** or ***) to get familiar with " + Session.getSession().getPrimaryModule().getDisplayName() + ". ", exampleLinks, VisualConstants.EXAMPLE_SESSION_ICON, c, this);
			}			
		}
	
		
		addLink("*** to continue working on previous sessions.", sessionLink, VisualConstants.OPEN_SESSION_LINK_ICON, c, this);

		
		// common links
		List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
		importLinks.add(importLink);
		importLinks.add(importFolderLink);
		importLinks.add(importURLLink);

		// module specific links
		if (!application.isStandalone) {
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
	
	
	private final static int MAX_ROW_CHARS = 42;

	public static void addLink(String description, JXHyperlink link, ImageIcon icon, GridBagConstraints c, JComponent component) {
	
		LinkUtil.addLink(description, link, icon, c, component, MAX_ROW_CHARS, Color.white);
	}

	
	public static void addLinks(String description, List<JXHyperlink> links, ImageIcon icon, GridBagConstraints c, JComponent component) {
		LinkUtil.addLinks(description, links, icon, c, component, MAX_ROW_CHARS, Color.white);
	}

}
