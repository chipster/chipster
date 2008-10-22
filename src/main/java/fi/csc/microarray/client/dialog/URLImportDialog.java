package fi.csc.microarray.client.dialog;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportUtils;

/**
 * Dialog for asking URL from user. It has an drop down menu which shows the 
 * most recently typed URLs. Checks if URL is not correctly formated and 
 * shows an error message if it is malformed. Does not check that the URL 
 * points to valid location.
 * 
 * @author mkoski
 *
 */
public class URLImportDialog extends JDialog implements ActionListener{
	
	private final Dimension BUTTON_SIZE = new Dimension(70, 25);
	private final Dimension COMBO_SIZE  = new Dimension(300, 25);
	private final Dimension LABEL_SIZE  = new Dimension(100, 25);
	
	private final int RECENT_URL_COUNT = 10;
	
	private static final String SKIP_TEXT = "Import directly if possible";
	private JCheckBox skipCheckBox;
	
	private JLabel label;
	private JButton okButton;
	private JButton cancelButton;
	private JComboBox URLComboBox;
	private JComboBox folderNameCombo;
	private static LinkedList<String> recentURLs;
	private SwingClientApplication client;
	private URL selectedURL;
	
	public URLImportDialog(SwingClientApplication client) {
		super(client.getMainFrame(), true);
		if(recentURLs == null){
			recentURLs = new LinkedList<String>();
		}
		
		this.client = client;
		this.setTitle("Import data from URL");
		this.setModal(true);
		
		label = new JLabel("Insert URL");
		label.setFont(label.getFont().deriveFont(Font.BOLD));
		label.setPreferredSize(LABEL_SIZE);
		
		folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true).toArray());
		folderNameCombo.setEditable(true);
		
		skipCheckBox = new JCheckBox(SKIP_TEXT);
		skipCheckBox.setSelected(true);
		
		okButton = new JButton("OK");
		okButton.setPreferredSize(BUTTON_SIZE);
		okButton.addActionListener(this);
		
		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);
		
		String[] comboValues = new String[recentURLs.size() + 1];
		comboValues[0] = "http://";
		for(int i = 0; i < recentURLs.size(); i++){
			comboValues[i+1] = recentURLs.get(i);
		}
		URLComboBox = new JComboBox(comboValues);
		URLComboBox.setBackground(Color.WHITE);
		URLComboBox.setPreferredSize(COMBO_SIZE);
		URLComboBox.setEditable(true);
				
		GridBagConstraints c = new GridBagConstraints();
		
		this.setLayout(new GridBagLayout());
		// Label
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(10, 10, 5, 10);
		c.gridx = 0; 
		c.gridy = 0;
		this.add(label, c);
		
		// Combobox
		c.insets.set(0,10,10,10);		
		c.gridy++;
		this.add(URLComboBox, c);
		
		c.insets.set(10,10,5,10);
		c.gridy++;
		this.add(new JLabel("Insert in folder"),c);
		
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets.set(0,10,10,10);
		c.gridy++;
		this.add(folderNameCombo,c);
		

		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.EAST;
		c.gridy++;
		this.add(skipCheckBox,c);
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridy++;
		JPanel keepButtonsRightPanel = new JPanel();
		keepButtonsRightPanel.add(okButton);
		keepButtonsRightPanel.add(cancelButton);
		this.add(keepButtonsRightPanel, c);
		
		this.pack();
		this.setLocationRelativeTo(null);
		URLComboBox.requestFocusInWindow();
		this.setVisible(true);
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == okButton){
			//Most  important functionality in the method where this class is created
			
			setSelectedURL();
			this.dispose();
		} else if (e.getSource() == cancelButton){
			selectedURL = null;
			this.dispose();
		}
	}
	
	/**
	 * Gets the recently selected URL
	 * @return
	 */
	public URL getSelectedURL(){
		return selectedURL;
	}
	
	public String getSelectedFolderName(){
		return folderNameCombo.getSelectedItem().toString();
	}
	
	public boolean isSkipSelected(){
		return skipCheckBox.isSelected();
	}

	/**
	 * Gets URL from combobox and sets it to selectedURL field. Displays an error message 
	 * dialog if url is. incorrectly formated. Adds the URL to recent typed URL list.
	 * 
	 * @return
	 */
	private void setSelectedURL(){
		URL url;
		try {
			url = new URL(URLComboBox.getSelectedItem().toString());
		} catch (MalformedURLException e1) {
			JOptionPane.showMessageDialog(client.getMainFrame(), "Malformed URL", "Malformed URL", JOptionPane.ERROR_MESSAGE);
			url = null;
		}
		
		if(url != null){
			selectedURL = url;
			addToRecentURLList(url.toString());
		} else {
			selectedURL = null;
		}
	}

	private void addToRecentURLList(String url){
		if(!recentURLs.contains(selectedURL.toString())){
			if(recentURLs.size() == RECENT_URL_COUNT){
				recentURLs.removeLast();
			} 
			recentURLs.addFirst(url);
		}
	}
}
