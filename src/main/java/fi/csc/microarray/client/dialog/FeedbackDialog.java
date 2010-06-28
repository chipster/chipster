package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;

/**
 * Dialog for user to give feedback, report an error etc.
 * 
 * FIXME: currently not used, should be shown in error dialogs
 * after clicking "Send Report".
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
        this.setTitle("Send Report");
        this.setModal(true);
        
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
        
        // Show
        this.pack();
        this.setResizable(false);
        this.setLocationRelativeTo(application.getMainFrame());
        this.setVisible(true);
    }

    @Override
    public void actionPerformed(ActionEvent event) {
        // check mandatory fields
        if (detailArea.getText().length() == 0) {
            application.showDialog("Please, provide some details.", Severity.INFO, true);
        }
        
        // TODO send feedback message using manager topic
    }
}
