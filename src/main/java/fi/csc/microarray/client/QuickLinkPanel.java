package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
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

public class QuickLinkPanel extends JPanel implements ActionListener{
	
	private SwingClientApplication application;

	private JXHyperlink sessionLink;
	private JXHyperlink importLink;
	private JXHyperlink emptyLink;
	private JXHyperlink exampleLink;
	private JXHyperlink importFolderLink;
	private JXHyperlink importURLLink;

	private static final String LINK_WORD = "***";

	public QuickLinkPanel() {
		super(new GridBagLayout());
		
		application = (SwingClientApplication)Session.getSession().getApplication();
		
		this.setBackground(Color.white);
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;		
		c.insets.set(5, 10, 5, 10);
		
		c.anchor = GridBagConstraints.NORTHWEST;
		c.gridwidth = 2;
		
		JLabel title = new JLabel("Getting started");
		title.setFont(title.getFont().deriveFont(
				(float) (title.getFont().getSize()*1.2)));
		title.setFont(title.getFont().deriveFont(Font.BOLD));
		
		this.add(title, c);
		c.gridwidth = 1;		
		c.insets.set(0, 10, 0, 0);
				
		addLink("*** to continue working on previous sessions.", getSessionLink(),
				VisualConstants.OPEN_SESSION_LINK_ICON,	c);

		List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
		importLinks.add(getImportLink());
		importLinks.add(getImportFolderLink());
		importLinks.add(getImportURLLink());
		
		addLink("Import new data to Chipster: \n    *** \n    *** \n    ***", importLinks,  
				VisualConstants.IMPORT_LINK_ICON, c);
								
		addLink("*** to learn more and play around.", getExampleLink(), 
				VisualConstants.EXAMPLE_SESSION_ICON, c);
		
		//Panels to take rest of space
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
		
		
		this.setMinimumSize(new Dimension(0,0));
		this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));
	}
	
	private void addLink(String description, JXHyperlink link, ImageIcon icon, GridBagConstraints c){
		List<JXHyperlink> list = new LinkedList<JXHyperlink>();
		list.add(link);
		addLink(description, list, icon, c);
	}
	
	private void addLink(String description, List<JXHyperlink> links, ImageIcon icon, GridBagConstraints c){
		
		String[] words = description.split(" ");
		int rowChars = 0;
		final int MAX_ROW_CHARS = 40;
		Iterator<JXHyperlink> linkIterator = links.iterator();
		int rowCount = 0;
			
		c.gridx = 1;
		c.insets.top = 10;
		JPanel row = null;
						
		for (int i = 0; i < words.length; i++){
			if(row == null || rowChars + words[i].length() > MAX_ROW_CHARS || words[i].equals("\n")){
				
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
						
			if(words[i].equals(LINK_WORD)){
				
				JXHyperlink link = linkIterator.next();
				rowChars += link.getText().length() + 1;
				row.add(link);
				//row.add(new JLabel(link.getText()));
			} else if(!words[i].equals("\n")) {	
				JLabel text = new JLabel(" " + words[i]);				
				row.add(text);
				rowChars += words[i].length() + 1;
			}
		}
		
		c.gridy -= (rowCount - 1);
		c.gridheight = rowCount;
		c.gridx = 0;
		c.insets.top = 10;
		this.add(new JLabel(icon), c);
		c.gridy += (rowCount - 1);
		c.gridheight = 1;
	}

	private JXHyperlink getExampleLink() { 
		if (exampleLink == null){
			exampleLink = new JXHyperlink();
			exampleLink.setText("Open example session");
			exampleLink.addActionListener(this);
		}
		return exampleLink;
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
