package fi.csc.microarray.client.dialog;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.Authenticator.LoginCallback;
import fi.csc.microarray.constants.VisualConstants;

@SuppressWarnings("serial")
public class LoginDialog extends JDialog implements ActionListener, KeyListener {

	private static final String LOGIN_ACTION = "login";
	private static final String CANCEL_ACTION = "cancel";
	private JTextField usernameField;
	private JPasswordField passwordField;
	private LoginCallback loginCallback;
	private boolean previousAttemptFailed = false;

	public static void main(String[] args) {
		LoginDialog dialog = new LoginDialog(null, true);
		dialog.setVisible(true);
	}

	public LoginDialog(LoginCallback callback) {
		this(callback, false);
	}

	public LoginDialog(LoginCallback callback, boolean previousAttemptFailed) {
		super();
		this.setMaximumSize(new Dimension(100, 100));
		this.loginCallback = callback;
		this.previousAttemptFailed = previousAttemptFailed;
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		setAlwaysOnTop(true);
		setResizable(false);
		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(getBanner(), BorderLayout.NORTH);
		getContentPane().add(getControls(), BorderLayout.CENTER);

		// Sets Look N Feel
		SwingClientApplication.setPlastic3DLookAndFeel(this);

		setLocationRelativeTo(null); // center on screen
		pack();
		addKeyListener(this);
	}

	private Component getControls() {
		JPanel panel = new JPanel();
		panel.setBorder(new EmptyBorder(40, 40, 10, 10));
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		panel.setLayout(gridbag);

		c.anchor = GridBagConstraints.WEST;
		c.insets = new Insets(2, 2, 12, 10);
		c.gridx = 0;
		c.gridy = 0;

		 c.gridwidth = 2;
		 
		 String announcement = Session.getSession().getApplication().getAnnouncementText();
		 String msg = 
				 "<html><p>Please enter your Chipster username and password,<br>"
				 + "or use the username 'guest' and password 'guest'<br/>to "
				 + "have a look. Running tools is disabled for guests.</p>";
		 if (announcement != null) {
			 msg += announcement;
		 }
		 msg  += "</html>";
		 
		 JLabel infoLabel = new JLabel(msg);
		 panel.add(infoLabel, c);
		 c.gridy++;
		 c.gridwidth = 1;

		if (previousAttemptFailed) {
			c.gridwidth = 2;
			JLabel failedLabel = new JLabel("Login failed, please check your username and password");
			failedLabel.setFont(failedLabel.getFont().deriveFont(Font.ITALIC));
			panel.add(failedLabel, c);
			c.gridwidth = 1;
			c.gridy++;
		}

		c.insets = new Insets(2, 2, 2, 10);
		panel.add(new JLabel("Username"), c);
		c.insets = new Insets(2, 2, 2, 2);
		c.gridx++;
		usernameField = new JTextField(24);
		usernameField.addKeyListener(this);
		panel.add(usernameField, c);
		c.gridy++;
		c.gridx = 0;

		c.insets = new Insets(2, 2, 2, 10);
		panel.add(new JLabel("Password"), c);
		c.insets = new Insets(2, 2, 2, 2);
		c.gridx++;
		passwordField = new JPasswordField(24);
		passwordField.addKeyListener(this);
		panel.add(passwordField, c);

		c.gridy++;
		c.gridx = 1;
		c.anchor = GridBagConstraints.EAST;
		c.insets = new Insets(15, 2, 0, 2);

		JPanel buttonPanel = new JPanel();

		JButton loginButton = new JButton("Login");
		loginButton.setActionCommand(LOGIN_ACTION);
		loginButton.addActionListener(this);
		loginButton.addKeyListener(this);
		buttonPanel.add(loginButton);
		JButton cancelButton = new JButton("Cancel");
		cancelButton.setActionCommand(CANCEL_ACTION);
		cancelButton.addActionListener(this);
		buttonPanel.add(cancelButton);

		panel.add(buttonPanel, c);

		return panel;
	}

	private Component getBanner() {
		return new JLabel(VisualConstants.getIcon(VisualConstants.LOGIN_BANNER));
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand() == LOGIN_ACTION) {
			login();

		} else if (e.getActionCommand() == CANCEL_ACTION) {
			loginCallback.cancel();
			dispose();
			//System.exit(1);

		} else {
			throw new RuntimeException("unknown action command: "
					+ e.getActionCommand());
		}

	}

	private void login() {
		loginCallback.authenticate(usernameField.getText(), new String(
				passwordField.getPassword()));
		dispose();
	}

	public void keyTyped(KeyEvent e) {
		if (e.getKeyChar() == '\n') {
			login();
		}
	}

	public void keyPressed(KeyEvent e) {
		// ignore
	}

	public void keyReleased(KeyEvent e) {
		// ignore
	}
}
