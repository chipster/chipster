package fi.csc.microarray.client;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URL;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.MicroarrayException;

public class QuickLinkPanel extends JPanel implements ActionListener{
	
	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink emptyLink;
	private JXHyperlink exampleLink;

	private CardLayout cardLayout;
	private JPanel cardPrarent;

	public QuickLinkPanel(JPanel parent, CardLayout cardLayout) {
		super(new GridBagLayout());
		
		this.cardLayout = cardLayout;
		this.cardPrarent = parent;
		
		application = (SwingClientApplication)Session.getSession().getApplication();
		
		this.setBackground(Color.white);
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;						
		c.insets.top = 5;
		c.insets.left = 10;
		//c.insets.right = 10;
		c.insets.bottom = 5;
		c.anchor = GridBagConstraints.LINE_START;
		
		JLabel title = new JLabel("Getting started");
		title.setFont(title.getFont().deriveFont(
				(float) (title.getFont().getSize()*1.2)));
		title.setFont(title.getFont().deriveFont(Font.BOLD));
		
		this.add(title,c);
		
		c.anchor = GridBagConstraints.LINE_START;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;
		c.weighty = 0.0;
		c.insets.bottom = 0;
		
		addLink(getSessionLink(), 
				"to continue analysis from the the saved situation.", c);
		addLink(getImportLink(), 
				"to analyse one or more data files saved in chipster \n" +
				"or another programs.", c);
		addLink(getExampleLink(), 
				"to try out Chipster analyses and visualisations.", c);
		//addLink(getEmptyLink(), "", c);
		
		this.setMinimumSize(new Dimension(0,0));
		this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));
	}
	
	private void addLink(Component link, String description, GridBagConstraints c){
		c.gridy++;
		c.insets.left = 10;
		c.insets.top = 5;
		this.add(link, c);
		
		if(description != null){
			c.gridy++;
			c.insets.left = 20;
			c.insets.top = 0;

			JTextArea text = new JTextArea(description);
			text.setEditable(false);

			this.add(text, c);
		}
	}

	private Component getExampleLink() { 
		if (exampleLink == null){
			exampleLink = new JXHyperlink();
			exampleLink.setText("Open example session");
			exampleLink.addActionListener(this);
		}
		return exampleLink;
	}

	private Component getEmptyLink() {
		if (emptyLink == null){
			emptyLink = new JXHyperlink();
			emptyLink.setText("Start empty session");
			emptyLink.addActionListener(this);
		}
		return emptyLink;
	}

	private Component getImportLink() {
		if (importLink == null){
			importLink = new JXHyperlink();
			importLink.setText("Import files");
			importLink.addActionListener(this);
		}
		return importLink;
	}

	private Component getSessionLink() {
		if (sessionLink == null){
			sessionLink = new JXHyperlink();
			sessionLink.setText("Open session");
			sessionLink.addActionListener(this);
		}
		return sessionLink;
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == sessionLink){
			application.loadSession();
		} else if(e.getSource() == importLink){
			try {
				application.openFileImport();
			} catch (Exception ex) {
				application.reportException(ex);
			}
		} else if(e.getSource() == emptyLink){
			
		} else if (e.getSource() == exampleLink){
			try {
				URL url = new URL("http://chipster.csc.fi/examples/kidney.cs");
				application.loadSessionFrom(url);
				
			} catch (Exception ex) {
				application.reportException(ex);
			}
		}
		
		cardLayout.last(cardPrarent);
	}
}
