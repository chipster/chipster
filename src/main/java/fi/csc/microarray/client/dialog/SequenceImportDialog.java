package fi.csc.microarray.client.dialog;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;

/**
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class SequenceImportDialog extends JDialog implements CaretListener, ActionListener {
    
    private String TASK_NAME = "\"Internal\"/\"importseq\"";
    private String OUT_FILE = "outseq";
    private final Dimension BUTTON_SIZE = new Dimension(70, 25);
    
    private static Logger logger = Logger.getLogger(SequenceImportDialog.class);
    private SwingClientApplication client;
    
    private JLabel nameLabel;
    private JLabel textLabel;
    private JTextField nameField;
    private JTextArea textArea;
    private JButton okButton;
    private JButton cancelButton;
    private JComboBox folderNameCombo;
    private JComboBox dbNameCombo;
    
    private enum Databases {
        PDB("PDB", "pdb"),
        EMBL("EMBL", "embl"),
        EMBL_NEW("EMBL New", "embl"),
        UNIPROT_SWISS("UniProt / Swiss", "swiss"),
        UNIPROT_TREMBL("UniProt / TrEMBL", "uniprot");
        
        private String name;
        private String value;
        
        private Databases(String name, String value) {
            this.name = name;
            this.value = value;
        }
        
        public String getValue() {
            return value;
        }
        
        public String toString() {
            return name;
        }
    }
    
    public SequenceImportDialog(SwingClientApplication client) {
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
        nameField = new JTextField();
        nameField.setPreferredSize(new Dimension(150, 20));
        nameField.setText("data.txt");
        c.insets.set(0,10,10,10);       
        c.gridy++;
        this.add(nameField, c);
        
        // Folder to store the file
        folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true).toArray());
        folderNameCombo.setPreferredSize(new Dimension(150, 20));
        folderNameCombo.setEditable(true);
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Create in folder"), c);
        c.insets.set(0,10,10,10);
        c.gridy++;
        this.add(folderNameCombo,c);
        
        // Database
        dbNameCombo = new JComboBox(Databases.values());
        dbNameCombo.setPreferredSize(new Dimension(150, 20));
        dbNameCombo.setBackground(Color.WHITE);
        c.insets.set(10,10,5,10);
        c.gridy++;
        this.add(new JLabel("Database"), c);
        c.insets.set(0,10,10,10);
        c.gridy++;
        this.add(dbNameCombo,c);
        
        // Text label
        textLabel = new JLabel("Sequence identifiers");
        c.anchor = GridBagConstraints.WEST;
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(textLabel, c);
        
        // Text area
        textArea = new JTextArea();
        textArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        textArea.setBorder(BorderFactory.createEmptyBorder());
        textArea.addCaretListener(this);
        JScrollPane areaScrollPane = new JScrollPane(textArea);
        areaScrollPane.setVerticalScrollBarPolicy(
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        areaScrollPane.setHorizontalScrollBarPolicy(
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        areaScrollPane.setBorder(BorderFactory.createLineBorder(Color.GRAY));
        areaScrollPane.setPreferredSize(new Dimension(150, 60));
        c.insets.set(0,10,10,10);  
        c.gridy++;
        this.add(areaScrollPane, c);
        
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
        this.setLocationRelativeTo(client.getMainFrame());
        this.setVisible(true);
    }

    public void caretUpdate(CaretEvent e) {
        
    }

    public void actionPerformed(ActionEvent e) {
        if(e.getSource() == okButton) {
            String fileName = this.nameField.getText();
            String folderName = (String) (this.folderNameCombo.getSelectedItem());
            String db = ((Databases) this.dbNameCombo.getSelectedItem()).getValue();
            String ids = this.textArea.getText();
            
            int separator = fileName.indexOf(".");
            separator = separator == -1 ? fileName.length() : separator;
            
            // There can be multiple sequences
            String[] sequences = ids.split("\n|;|,");
            int index = 0;
            for (String id : sequences) {
                // Rename files if we have to import several
                String name = fileName;
                if (sequences.length > 1) {
                    name = fileName.substring(0, separator) + "." +
                           index + fileName.substring(separator, fileName.length());
                    index++;
                }
                
                // Run import
                runImport(name, folderName, db, id);
            }
            this.dispose();
        } else if (e.getSource() == cancelButton) {
            this.dispose();
        }
    }
    
    /**
     * Run import job.
     * 
     * @param fileName
     * @param folderName
     * @param db
     * @param id
     */
    private void runImport(String fileName, String folderName, String db, String id) {
        try {               
            // Create importing job
            logger.info("Importing sequence...");
            TaskExecutor taskExecutor = new TaskExecutor(client.getEndpoint(), client.getDataManager());
            final Task importSequence = taskExecutor.createTask(TASK_NAME, true);
            importSequence.addParameter("sequence", db + ":" + id);
            
            // Run the job (blocking while it is progressing)
            taskExecutor.execute(importSequence);
            
            // Create a dataset
            DataBean data = importSequence.getOutput(OUT_FILE);
            data.setName(fileName);
            data.setContentType(client.getDataManager().guessContentType(fileName));
            data.setOperation(new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data }));
            
            // Make it visible
            DataFolder folder = client.initializeFolderForImport(folderName);
            folder.addChild(data);
            client.getSelectionManager().selectSingle(data, this);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    }
}
