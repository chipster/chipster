package fi.csc.microarray.client.dialog;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.LinkedList;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
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
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class SequenceImportDialog extends JDialog implements CaretListener, ActionListener {
    
    private static final String OPERATION_ID = "import.sadl";
    private static final String OUT_FILE = "outseq";
    private final Dimension BUTTON_SIZE = new Dimension(70, 25);
    
    private static Logger logger = Logger.getLogger(SequenceImportDialog.class);
    private SwingClientApplication application;
    
    private JLabel nameLabel;
    private JLabel textLabel;
    private JTextField nameField;
    private JTextField beginField;
    private JTextField endField;
    private JTextArea textArea;
    private JButton okButton;
    private JButton cancelButton;
    private JComboBox folderNameCombo;
    private JComboBox dbNameCombo;
    private JCheckBox mergeCheckBox;
    
    private enum Databases {
        EMBL("EMBL", "embl"),
        EMBL_NEW("EMBL New", "emblnew"),
        PDB("PDB", "pdb_seq"),
        UNIPROT_SWISS("UniProt / SwissProt", "swiss"),
        UNIPROT_TREMBL("UniProt / TrEMBL", "trembl");
        
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
    
    public SequenceImportDialog(SwingClientApplication application) {
        super(application.getMainFrame(), true);

        this.application = application;
        this.setTitle("Import sequence");
        this.setModal(true);
        
        // Layout
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.WEST;
        c.gridx = 0; 
        c.gridy = 0;
        this.setLayout(new GridBagLayout());
        
        // Name label
        nameLabel = new JLabel("Filename");
        c.anchor = GridBagConstraints.WEST;
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(nameLabel, c);
        
        // Name field
        nameField = new JTextField();
        nameField.setPreferredSize(new Dimension(150, 20));
        nameField.setText("data.txt");
        c.insets.set(0, 10, 10, 10);       
        c.gridy++;
        this.add(nameField, c);
        
        // Folder to store the file
        folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true).toArray());
        folderNameCombo.setPreferredSize(new Dimension(150, 20));
        folderNameCombo.setEditable(true);
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(new JLabel("Create in folder"), c);
        c.insets.set(0, 10, 10, 10);
        c.gridy++;
        this.add(folderNameCombo, c);
        
        // Database
        dbNameCombo = new JComboBox(Databases.values());
        dbNameCombo.setPreferredSize(new Dimension(150, 20));
        dbNameCombo.setBackground(Color.WHITE);
        dbNameCombo.setSelectedItem(Databases.UNIPROT_SWISS);
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
        areaScrollPane.setPreferredSize(new Dimension(400, 90));
        c.insets.set(0, 10, 10, 10);  
        c.gridy++;
        this.add(areaScrollPane, c);
        
        // Range fields
        beginField = new JTextField(3);
        endField = new JTextField(3);
        JPanel rangePanel = new JPanel();
        rangePanel.setLayout(new FlowLayout(FlowLayout.LEFT));
        rangePanel.add(new JLabel("Start"));
        rangePanel.add(beginField);
        rangePanel.add(new JLabel("End"));
        rangePanel.add(endField);
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(rangePanel, c);
        
        // Checkbox for merging into one
        mergeCheckBox = new JCheckBox("Merge into one dataset");
        mergeCheckBox.setSelected(true);
        c.insets.set(10, 10, 5, 10);
        c.gridy++;
        this.add(mergeCheckBox, c);
        
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
        
        // Default focus
        textArea.requestFocusInWindow();
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
            LinkedList<DataBean> datasets = new LinkedList<DataBean>();
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
                datasets.add(runImport(name, folderName, db, id));
            }
            
            // Create datasets
            if (mergeCheckBox.isSelected() && datasets.size() > 0) {
                // Merge all imported datasets into one
                try {
                    StringBuilder sb = new StringBuilder();
                    for (DataBean dataset : datasets) {
                        sb.append(new String(dataset.getContents()));
                    }
                    String content = sb.toString();
                    
                    // Create single dataset for all input data
                    ByteArrayInputStream stream = new ByteArrayInputStream(content.getBytes());
                    DataBean data = application.getDataManager().createDataBean(fileName, stream);
                    data.setContentType(application.getDataManager().guessContentType(fileName));
                    data.setOperation(new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data }));
                    
                    // Make it visible
                    DataFolder folder = application.initializeFolderForImport(folderName);
                    folder.addChild(data);
                    application.getSelectionManager().selectSingle(data, this);
                } catch (MicroarrayException e1) {
					// FIXME proper error handling 
                	e1.printStackTrace();
                } catch (IOException ioe) {
					// FIXME proper error handling 
					ioe.printStackTrace();
				}
            } else {
                // Make separate datasets visible
                for (DataBean dataset : datasets) {
                    DataFolder folder = application.initializeFolderForImport(folderName);
                    folder.addChild(dataset);
                    application.getSelectionManager().selectSingle(dataset, this);
                }
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
    private DataBean runImport(String fileName, String folderName, String db, String id) {
        DataBean data = null;
        try {               
            // Create importing job
            logger.info("Importing sequence...");
            
            Operation operation = new Operation(application.getOperationDefinition(OPERATION_ID), new DataBean[] {});
            operation.setParameter("sequence", db + ":" + id);
            operation.setParameter("sbegin", beginField.getText());
            operation.setParameter("send", endField.getText());
            
            // Run the job (blocking while it is progressing)
            application.executeOperation(operation);

            // should not be done like this
//            // Create a dataset or prepare for merging them later
//            data = importSequence.getOutput(OUT_FILE);
//            data.setName(fileName);
//            data.setContentType(client.getDataManager().guessContentType(fileName));
//            data.setOperation(new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data }));
        } catch (Exception exc) {
            // FIXME proper error handling
        	exc.printStackTrace();
        }
        return data;
    }
}
