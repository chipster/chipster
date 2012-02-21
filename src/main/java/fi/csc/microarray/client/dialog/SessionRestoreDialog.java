package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Date;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.SwingClientApplication;

/**
 * A modal dialog which allows user to restore a failed session. 
 * 
 */
public class SessionRestoreDialog extends JDialog implements ActionListener {

	private final Dimension BUTTON_SIZE = new Dimension(75, 25);
	
	private JButton restoreButton;
	private JButton discardButton;
	
	private ClientApplication client;
	private File sessionFile;
	
	/**
	 * Creates a new data item rename dialog and shows it.
	 * 
	 * @param owner The owning frame component of the dialog to be created.
	 * 				In this case, it would be the GUI mainframe class.
	 * @param sessionFile	session file to restore
	 */
	public SessionRestoreDialog(SwingClientApplication client, File sessionFile) {
		super((JFrame)null/*client.getMainFrame()*/, "Restore session");
		
		this.client = client;
		this.sessionFile = sessionFile;
		
		restoreButton = new JButton("Restore");
		restoreButton.setPreferredSize(BUTTON_SIZE);
		discardButton = new JButton("Discard");
		discardButton.setPreferredSize(BUTTON_SIZE);
		
		JPanel contentPane = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.WEST;
		c.insets = new Insets(10, 20, 5, 10);
		c.gridx = 0; c.gridy = 0;
		c.gridwidth = 2;
		JLabel dataSetNameLabel = new JLabel("<html>It seems that Chipster crashed " + new Date(sessionFile.lastModified()) + ".<br>Do you want to restore the data you were working on?</html>");
		contentPane.add(dataSetNameLabel, c);
		c.gridy = 1;
		c.gridwidth = 1;
		contentPane.add(restoreButton, c);
		c.gridx = 1;
		contentPane.add(discardButton, c);
		this.setContentPane(contentPane);
		
		restoreButton.addActionListener(this);
		discardButton.addActionListener(this);
		
		this.pack();
		this.setLocationRelativeTo(null);  // centered on screen
		this.setVisible(true);
	}
	
	/**
	 * A method defined by the Actionlistener interface that allows this
	 * dialog to react accordingly when its buttons are clicked.
	 * 
	 * @param e An event heralding a click on one of this dialog's buttons.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		
		if (source == restoreButton) {
			client.restoreSessionFrom(sessionFile); // client will clear dead dirs after the restore 
			
		} else {
			client.clearDeadTempDirectories();
		}
		
		this.dispose();
	}
}