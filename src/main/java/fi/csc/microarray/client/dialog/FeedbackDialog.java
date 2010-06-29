package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.jms.JMSException;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.FeedbackMessage;

/**
 * Dialog for user to give feedback, report an error etc.
 * 
 * @author naktinis
 *
 */
public class FeedbackDialog extends JDialog implements ActionListener {
    
    private SwingClientApplication application;
    private final Dimension BUTTON_SIZE = new Dimension(70, 25);
    
    private JButton okButton;
    private JButton cancelButton;
    private JTextArea detailArea;
    private JTextField emailField;
    
    public FeedbackDialog(SwingClientApplication application) {
        super(application.getMainFrame(), true);

        this.application = application;
        this.setTitle("Send a report");
        
        // Layout
        this.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.WEST;
        c.gridx = 0; 
        c.gridy = 0;
        
        // Text are for entering details
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Details"), c);
        detailArea = new JTextArea();
        detailArea.setPreferredSize(new Dimension(300, 150));
        c.insets.set(0, 10, 10, 10);  
        c.gridy++;
        this.add(detailArea, c);
        
        // Email
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Email"), c);
        emailField = new JTextField();
        emailField.setPreferredSize(new Dimension(300, 20));
        c.insets.set(0, 10, 10, 10);  
        c.gridy++;
        this.add(emailField, c);
        
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
                        "Please, provide some details.",
                        "Info", JOptionPane.INFORMATION_MESSAGE);
            }
            
            // send feedback message using manager topic
            try {
                FeedbackMessage message = new FeedbackMessage(detailArea.getText(), emailField.getText());
                MessagingTopic requestTopic = application.getEndpoint().createTopic(Topics.Name.FEEDBACK_TOPIC,
                        AccessMode.WRITE);
                // TODO store interesting information in the file broker
                // create a temporary session file
                // File.createTempFile(arg0, arg1)
                // getDataManager().saveSnapshot(file, application);
                // save it with the file broker
                // this.fileBroker = new FileBrokerClient(endpoint.createTopic(Topics.Name.URL_TOPIC, AccessMode.WRITE));
                // fileBroker.addFile(bean.getContentByteStream(), progressListener);
                // TODO delete temp file
                requestTopic.sendMessage(message);
                this.dispose();
            } catch (JMSException e) {
                application.reportException(e);
            }
            
        } else {
            // close
            this.dispose();
        }
    }
}
