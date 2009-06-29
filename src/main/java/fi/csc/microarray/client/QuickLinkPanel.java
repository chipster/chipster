package fi.csc.microarray.client;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;

import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

public class QuickLinkPanel extends JPanel implements ActionListener{
	
	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink emptyLink;
	private JXHyperlink exampleLink;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;

	private CardLayout cardLayout;
	private JPanel cardParent;


	public QuickLinkPanel(JPanel parent, CardLayout cardLayout) {
		super(new GridBagLayout());
		
		this.cardLayout = cardLayout;
		this.cardParent = parent;
		
		application = (SwingClientApplication)Session.getSession().getApplication();
		
		this.setBackground(Color.white);
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;						
		c.insets.top = 5;
		c.insets.left = 10;
		//c.insets.right = 10;
		c.insets.bottom = 5;
		c.anchor = GridBagConstraints.NORTHWEST;
		
		JLabel title = new JLabel("Getting started");
		title.setFont(title.getFont().deriveFont(
				(float) (title.getFont().getSize()*1.2)));
		title.setFont(title.getFont().deriveFont(Font.BOLD));
		
		this.add(title,c);
		
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;
		c.weighty = 0.0;
		c.insets.bottom = 0;
		
		addLink(getSessionLink(), 
				"to continue working on previous sessions.", c);
		
		JXHyperlink emptyLink = new JXHyperlink();
		emptyLink.setText("");
		addLink(emptyLink, 
				"Import new data to Chipster:", c);
		
		c.gridy++;
		c.insets.top = 5;
		c.insets.left += 10;
		
		this.add(getImportLink(), c);
		c.gridy++;
		c.insets.top = 0;
		this.add(getImportFolderLink(), c);
		c.gridy++;
		this.add(getImportURLLink(), c);
		
		c.insets.left -= 10;
						
		addLink(getExampleLink(), 
				"to learn more and play around.", c);
		//addLink(getEmptyLink(), "", c);
		
		c.weighty = 1.0;
		c.fill = GridBagConstraints.BOTH;
		JPanel emptyPanel = new JPanel();
		emptyPanel.setBackground(Color.white);
		this.add(emptyPanel, c);
		
		this.setMinimumSize(new Dimension(0,0));
		this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));
	}
	
	private void addLink(JXHyperlink link, String description, GridBagConstraints c){
		
		String[] words = description.split(" ");
		int rowChars = link.getText().length() + 1;
		final int MAX_ROW_CHARS = 40;
			
		c.gridy++;
		c.insets.left = 10;
		c.insets.top = 15;
		c.insets.bottom = 0;
		
		JPanel row = null;
						
		for (int i = -1; i < words.length; i++){
			if(i == -1 || rowChars + words[i].length() > MAX_ROW_CHARS){
				
				FlowLayout flow = new FlowLayout(FlowLayout.LEADING);
				flow.setVgap(0);
				flow.setHgap(0);				
				row = new JPanel(flow);
				row.setBackground(Color.white);

				c.gridy++;
				this.add(row, c);
				c.insets.top = 0; // After first row
				
				if(i == -1){
					row.add(link);
				} else {
					rowChars = 0;					
				}
			} 
			
			if(i != -1){
				JLabel text = new JLabel(" " + words[i]);				
				row.add(text);
				rowChars += words[i].length();
			}
		}
	}

	private JXHyperlink getExampleLink() { 
		if (exampleLink == null){
			exampleLink = new JXHyperlink();
			exampleLink.setText("Open example session");
			exampleLink.addActionListener(this);
		}
		return exampleLink;
	}

	private JXHyperlink getEmptyLink() {
		if (emptyLink == null){
			emptyLink = new JXHyperlink();
			emptyLink.setText("Start empty session");
			emptyLink.addActionListener(this);
		}
		return emptyLink;
	}

	private JXHyperlink getImportLink() {
		if (importLink == null){
			importLink = new JXHyperlink();
			importLink.setText("Import files");
			importLink.addActionListener(this);
		}
		return importLink;
	}
	
	private JXHyperlink getImportFolderLink() {
		if (importFolderLink == null){
			importFolderLink = new JXHyperlink();
			importFolderLink.setText("Import folder");
			importFolderLink.addActionListener(this);
		}
		return importFolderLink;
	}
	
	private JXHyperlink getImportURLLink() {
		if (importURLLink == null){
			importURLLink = new JXHyperlink();
			importURLLink.setText("Import from URL");
			importURLLink.addActionListener(this);
		}
		return importURLLink;
	}

	private JXHyperlink getSessionLink() {
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
		} else if(e.getSource() == importURLLink){
			try {
				application.openURLImport();
			} catch (Exception ex) {
				application.reportException(ex);
			}
		} else if(e.getSource() == importFolderLink){
			application.openDirectoryImportDialog();
		} else if(e.getSource() == emptyLink){
			
		} else if (e.getSource() == exampleLink){
			try {
				URL url = new URL("http://chipster.csc.fi/examples/kidney.cs");
				application.loadSessionFrom(url);
				
			} catch (Exception ex) {
				application.reportException(ex);
			}
		}
	}
}
