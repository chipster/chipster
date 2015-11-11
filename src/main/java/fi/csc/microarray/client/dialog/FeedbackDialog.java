package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.jms.JMSException;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.messaging.message.FeedbackMessage;
import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.IOUtils;

/**
 * Dialog for user to give feedback, report an error etc.
 * 
 * @author naktinis
 *
 */
public class FeedbackDialog extends JDialog implements ActionListener {
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(FeedbackDialog.class);
	
    private SwingClientApplication application;
    private final Dimension BUTTON_SIZE = new Dimension(70, 25);
    
    private JButton okButton;
    private JButton cancelButton;
    private JTextArea detailArea;
    private JTextField emailField;
    private JCheckBox attachSessionBox;
    private JCheckBox attachLogsBox;
	private String errorMessage;

    public FeedbackDialog(SwingClientApplication application, String errorMessage) {
    	this(application, errorMessage, false);
    }
    	
    public FeedbackDialog(SwingClientApplication application, String errorMessage, boolean sendSessionByDefault) {
        super(application.getMainFrame(), true);

        this.application = application;
        this.setTitle("Contact support");
        this.errorMessage = errorMessage;
        
        // Layout
        this.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.WEST;
        c.gridx = 0; 
        c.gridy = 0;
        
        // Text are for entering details
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Message"), c);
        detailArea = new JTextArea();
        detailArea.setLineWrap(true);
        detailArea.setWrapStyleWord(true);
        detailArea.setPreferredSize(new Dimension(300, 150));
        c.insets.set(0, 10, 10, 10);  
        c.gridy++;
        this.add(detailArea, c);
        
        // Email
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Your email"), c);
        emailField = new JTextField();
        emailField.setPreferredSize(new Dimension(300, 20));
        c.insets.set(0, 10, 10, 10);  
        c.gridy++;
        this.add(emailField, c);
        
        // Checkbox for attaching user data
        c.insets.set(10,10,5,10);
        c.gridy++;
        attachSessionBox = new JCheckBox("Attach data and workflow information");
        attachSessionBox.setSelected(sendSessionByDefault);
        this.add(attachSessionBox, c);
        
        // Checkbox for client logs
        c.insets.set(10,10,5,10);
        c.gridy++;
        attachLogsBox = new JCheckBox("Attach log files");
        attachLogsBox.setSelected(true);
        this.add(attachLogsBox, c);
        
        // OK button
        okButton = new JButton("OK");
        okButton.setPreferredSize(BUTTON_SIZE);
        okButton.addActionListener(this);
        
        // Cancel button
        cancelButton = new JButton("Cancel");
        cancelButton.setPreferredSize(BUTTON_SIZE);
        cancelButton.addActionListener(this);
        
        // Buttons pannel
        JPanel buttonsPanel = new JPanel();
        buttonsPanel.add(okButton);
        buttonsPanel.add(cancelButton);
        c.fill = GridBagConstraints.HORIZONTAL;
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(buttonsPanel, c);
    }
    
    public void showDialog() {
        this.setLocationRelativeTo(application.getMainFrame());
        this.setResizable(false);
        this.pack();
        this.setModal(true);
        this.setVisible(true);
    }

    @Override
    public void actionPerformed(ActionEvent event) {
    	
        if (event.getSource() == okButton) {
            // check mandatory fields
            if (detailArea.getText().length() == 0) {
                JOptionPane.showMessageDialog(this.getParent(),
                        "Please, provide some details in the feedback area.",
                        "Info", JOptionPane.INFORMATION_MESSAGE);
                return;
            }

            this.dispose();
            
    		application.runBlockingTask("sending feedback, with large session this may take a while", new Runnable() {
    			public void run() {						
   
    	            // send feedback message using manager topic
    	            try {              
    	            	sendFeedback();
    	                
    	            } catch (Exception e) {
    	                application.reportException(e);
    	            }    			    			    			    
    			}
    		});
        } else {
            // close
            this.dispose();
        }
    }
    
    private void sendFeedback() throws Exception, JMSException, FileBrokerException, IOException, NotEnoughDiskSpaceException, FileNotFoundException {
    	ServiceAccessor serviceAccessor = Session.getSession().getServiceAccessor();
    	FileBrokerClient fileBrokerClient = serviceAccessor.getFileBrokerClient();

    	String sessionUrl;
    	if (attachSessionBox.isSelected()) {

    		String sessionId = application.getSessionManager().saveFeedbackSession();
    		sessionUrl = fileBrokerClient.getExternalURL(sessionId);

    	} else {
    		// user does not want to give data
    		sessionUrl = "";
    	}

    	// prepare feedback message
    	String messageText = detailArea.getText() + "\n\nError message:\n" + errorMessage;
    	FeedbackMessage message = new FeedbackMessage(messageText,
    			emailField.getText(), sessionUrl);

    	// attach log files if client allows
    	if (attachLogsBox.isSelected()) {
    		File logDir = DirectoryLayout.getInstance().getLogsDir();
    		for (File logFile : logDir.listFiles()) {
    			//create a copy of log, because otherwise it may grow during the upload and server side size check fails
    			File copyOfLog = File.createTempFile("chipster-temp-file", "");
    			copyOfLog.deleteOnExit();
    			IOUtils.copy(logFile, copyOfLog);

    			// save it with the file broker    	                    	
    			String dataId = CryptoKey.generateRandom();
    			fileBrokerClient.addFile(dataId, FileBrokerArea.CACHE, new FileInputStream(copyOfLog), copyOfLog.length(), null);

    			String url = fileBrokerClient.getExternalURL(dataId);
    			message.addLog(logFile.getName(), url);
    			copyOfLog.delete();
    		}
    	}

    	// send feedback message to manager
    	serviceAccessor.sendFeedbackMessage(message);
    }
}
