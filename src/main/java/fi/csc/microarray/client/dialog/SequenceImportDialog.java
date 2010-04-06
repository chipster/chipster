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

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataBean;

/**
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class SequenceImportDialog extends JDialog implements CaretListener, ActionListener {
    
    private String TASK_NAME = "\"Edit\"/\"seqret\"";
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
        PDB("PDB"),
        EMBL("EMBL"),
        EMBL_NEW("EMBL New"),
        UNIPROT_SWISS("UniProt / Swiss"),
        UNIPROT_TREMBL("UniProt / TrEMBL");
        
        private String name;
        
        private Databases(String name) {
            this.name = name;
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
        nameField.addCaretListener(this);
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
        if(e.getSource() == okButton){
            String fileName = this.nameField.getText();
            String fileContent = this.textArea.getText();
            String folderName = (String) (this.folderNameCombo.getSelectedItem());
            try {                
                // TODO
                
                // Create importing job
                logger.info("Importing sequence...");
                TaskExecutor taskExecutor = new TaskExecutor(client.getEndpoint(), client.getDataManager());
                final Task importSequence = taskExecutor.createTask(TASK_NAME, true);
                importSequence.addParameter("sequence", "swiss:CASA1_HUMAN");
                importSequence.addParameter("feature", "Y");
                importSequence.addParameter("firstonly", "N");
                //importSequence.addInput(name, input)("sequence", null);
                
                // Run the job (blocking while it is progressing)
                taskExecutor.execute(importSequence);
                
                // Parse metadata
                DataBean metadataBean = importSequence.getOutput(OUT_FILE);
                //this.metadata = new String(metadataBean.getContents());
                //manager.delete(metadataBean); // don't leave the bean hanging around
                //logger.debug("got metadata: " + this.metadata.substring(0, 50) + "...");
                //List<SADLDescription> descriptions = new ChipsterSADLParser().parseMultiple(this.metadata);
                //this.parsedCategories = new OperationGenerator().generate(descriptions).values();
            } catch (Exception exc) {
                exc.printStackTrace();
            }

            this.dispose();
        } else if (e.getSource() == cancelButton) {
            this.dispose();
        }
    }
}
