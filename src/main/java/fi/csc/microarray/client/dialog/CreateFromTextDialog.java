package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;

/**
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class CreateFromTextDialog extends JDialog implements CaretListener, ActionListener {
    
    private final Dimension BUTTON_SIZE = new Dimension(70, 25);
    
    private SwingClientApplication client;
    
    private JLabel nameLabel;
    private JLabel textLabel;
    private JTextField nameField;
    private JTextArea textArea;
    private JButton okButton;
    private JButton cancelButton;
    private JComboBox folderNameCombo;
    
    public CreateFromTextDialog(SwingClientApplication client) {
        super(client.getMainFrame(), true);

        this.client = client;
        this.setTitle("Create dataset from text");
        this.setModal(true);
        
        // Layout
        GridBagConstraints c = new GridBagConstraints();
        this.setLayout(new GridBagLayout());
        
        // Name label
        nameLabel = new JLabel("Filename");
        c.anchor = GridBagConstraints.WEST;
        c.insets.set(10, 10, 5, 10);
        c.gridx = 0; 
        c.gridy = 0;
        this.add(nameLabel, c);
        
        // Name field
        nameField = new JTextField(30);
        nameField.setText("data.txt");
        nameField.addCaretListener(this);
        c.insets.set(0,10,10,10);       
        c.gridy++;
        this.add(nameField, c);
        
        // Folder to store the file
        folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true).toArray());
        folderNameCombo.setEditable(true);
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Create in folder"), c);
        c.insets.set(0,10,10,10);
        c.gridy++;
        this.add(folderNameCombo,c);
        
        // Text label
        textLabel = new JLabel("Text");
        c.anchor = GridBagConstraints.WEST;
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(textLabel, c);
        
        // Text area
        textArea = new JTextArea(18, 43);
        textArea.setLineWrap(true);
        textArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        textArea.addCaretListener(this);
        c.insets.set(0,10,10,10);  
        c.gridy++;
        this.add(textArea, c);
        
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
        this.setLocationRelativeTo(client.getMainFrame());
        this.setVisible(true);
    }

    public void caretUpdate(CaretEvent e) {
        
    }

    public void actionPerformed(ActionEvent e) {
        if(e.getSource() == okButton){
            String fileName = this.nameField.getText();
            String fileContent = this.textArea.getText();
            String folderName = (String) (this.folderNameCombo.getSelectedItem());
            try {
                // Create temporary
                ByteArrayInputStream stream = new ByteArrayInputStream(fileContent.getBytes());

                // Create dataset
                DataBean data = client.getDataManager().createDataBean(fileName, stream);
                data.setContentType(client.getDataManager().guessContentType(fileName));
                data.setOperation(new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data }));
                
                // Make it visible
                DataFolder folder = client.initializeFolderForImport(folderName);
                folder.addChild(data);
                
                // Select
                client.getSelectionManager().selectSingle(data, this);
            } catch (Exception exc) {
                exc.printStackTrace();
            }

            this.dispose();
        } else if (e.getSource() == cancelButton){
            this.dispose();
        }
    }
}
