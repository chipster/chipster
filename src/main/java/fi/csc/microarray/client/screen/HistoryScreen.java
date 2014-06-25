package fi.csc.microarray.client.screen;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.util.GeneralFileFilter;

/**
 * A screen for viewing the history (longer version) of a dataset as text.
 * Allows the user to copypaste it to another application or save it as
 * a text file.
 * 
 * @author Janne Käki
 *
 */
public class HistoryScreen extends ScreenBase
                           implements ActionListener, CaretListener {
	private static final Dimension BUTTON_SIZE = new Dimension(100, 25);
	
	private static final ClientApplication application = Session.getSession().getApplication();
	
	private JFrame frame;
	private JTextArea textArea;
	private Map<String, JCheckBox> checkBoxes= new HashMap<String, JCheckBox>();
	private JButton saveButton;
	private JButton closeButton;
	private DataBean data;	
	private JFileChooser chooser = null;
	
	/**
	 * Initializes a new history screen.
	 */
	public HistoryScreen() {
		frame = new JFrame("History");
		JPanel contentPane = new JPanel(new GridBagLayout());
		
		checkBoxes.put("title", new JCheckBox("Step title"));
		checkBoxes.put("name", new JCheckBox("Dataset name"));
		checkBoxes.put("date", new JCheckBox("Creation date"));
		checkBoxes.put("oper", new JCheckBox("Applied analysis tool"));
		checkBoxes.put("param", new JCheckBox("Parameters"));
		checkBoxes.put("notes", new JCheckBox("User notes"));
		checkBoxes.put("code", new JCheckBox("Source code"));

		for (JCheckBox box : checkBoxes.values()) {
			box.setSelected(true);
			box.addActionListener(this);
		}
		
		checkBoxes.get("notes").setSelected(false);
		checkBoxes.get("date").setSelected(false);
		checkBoxes.get("code").setSelected(false);
		checkBoxes.get("param").setEnabled(checkBoxes.get("oper").isSelected());
        checkBoxes.get("code").setEnabled(checkBoxes.get("oper").isSelected() && !application.isStandalone());
		
		saveButton = new JButton("Save...");
		saveButton.setPreferredSize(BUTTON_SIZE);
		saveButton.addActionListener(this);
		
		closeButton = new JButton("Close");
		closeButton.setPreferredSize(BUTTON_SIZE);
		closeButton.addActionListener(this);
		
		textArea = new JTextArea();
		textArea.setText(getHistoryText());
		textArea.setMargin(new Insets(2, 2, 2, 2));
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		textArea.addCaretListener(this);
		
		JScrollPane textAreaScroller = new JScrollPane(textArea);
		textAreaScroller.setPreferredSize(new Dimension(500, 250));	
		
		JLabel topLabel = new JLabel("Show for Datasets:");
		
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(10, 2, 1, 2);
		c.gridx = 0; c.gridy = 0;
		contentPane.add(topLabel, c);
		c.insets.set(1, 2, 0, 10);
		c.gridwidth = 1;
		c.gridy++;
		contentPane.add(checkBoxes.get("title"), c);
		c.gridy++;
		contentPane.add(checkBoxes.get("name"), c);
		c.gridy++;
		contentPane.add(checkBoxes.get("date"), c);
		c.gridx++; c.gridy = 1;
		contentPane.add(checkBoxes.get("oper"), c);
		c.gridy++;
		//c.insets.set(1, 20, 0, 2);
		contentPane.add(checkBoxes.get("param"), c);
        c.gridy++;
      	contentPane.add(checkBoxes.get("code"), c);
        c.gridx++; c.gridy = 1;
        c.insets.set(1, 2, 0, 2);
        contentPane.add(checkBoxes.get("notes"), c);
		c.insets.set(10, 2, 1, 2);
		c.gridx = 0; c.gridy += 3;
		c.gridwidth = 3;
		c.weightx = 1;
		c.weighty = 1;
		c.fill = GridBagConstraints.BOTH;		
		contentPane.add(textAreaScroller, c);
		c.insets.set(2, 2, 2, 2);
		c.gridx = 0; c.gridy++;
		c.gridwidth = 1;
		c.weightx = 0;
		c.weighty = 0;
		c.fill = GridBagConstraints.NONE;
		contentPane.add(saveButton, c);
		c.gridx++;
		c.anchor = GridBagConstraints.EAST;
		contentPane.add(closeButton, c);
		
		frame.setContentPane(contentPane);
		frame.pack();
		frame.setResizable(true);
		frame.setLocationRelativeTo(null);  // centered on screen
	}
	
	/**
	 * Sets the data whose history is to be viewed.
	 * 
	 * @param data
	 */
	public void setData(DataBean data) {
		this.data = data;		
		refreshText();
	}
	
	private void refreshText() {
		textArea.setText(getHistoryText());
	}
	
	public String getHistoryText() {	
		String history = application.getHistoryText(data, 
				checkBoxes.get("title").isSelected(),
				checkBoxes.get("name").isSelected(),
				checkBoxes.get("date").isSelected(),
				checkBoxes.get("oper").isSelected(),
				checkBoxes.get("code").isSelected(),
				checkBoxes.get("notes").isSelected(),
				checkBoxes.get("param").isSelected());
		saveButton.setEnabled(true);
		return history;
	}
	
	/**
	 * Asks the user for a filename and saves the history text in that file.
	 * If one with that name already exists, the program will ask for
	 * confirmation. The saved text is equal to the one shown in the textarea.
	 * 
	 * @return True if save succeeded, false if the user aborted it or
	 * 		   if the writing of the file otherwise failed.
	 */
	private boolean exportToText() {
		
		if (chooser == null) {
			chooser = ImportUtils.getFixedFileChooser();
			GeneralFileFilter filter = new GeneralFileFilter("Text Files", new String[] { "txt" } );
			chooser.setFileFilter(filter);
		}

		int ret = chooser.showSaveDialog(null);
		if (ret == JFileChooser.APPROVE_OPTION ) {
			
			File selectedFile = chooser.getSelectedFile();
			if (selectedFile.exists()) {
				if (JOptionPane.showConfirmDialog(frame,
						"Really overwrite " + selectedFile.getName() + "?",
						"Save",
						JOptionPane.YES_NO_OPTION) != JOptionPane.YES_OPTION) {
					return false;  // overwrite was not confirmed
				}
			}
			
			try {
				PrintWriter writer = new PrintWriter(selectedFile);
				String text = textArea.getText();
				String[] lines = text.split("\n");
				for (String line : lines) {
					writer.print(line);
					writer.println();
				}
				writer.flush();
				writer.close();
				return true;
			} catch (IOException e) {
				application.reportException(e);
				return false;
			}
		}
		return false;
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source instanceof JCheckBox) {
			if (source == checkBoxes.get("oper")) {
				checkBoxes.get("param").setEnabled(checkBoxes.get("oper").isSelected());
                checkBoxes.get("code").setEnabled(checkBoxes.get("oper").isSelected() && !application.isStandalone());
			}
			textArea.setText(getHistoryText());
		} else if (source == saveButton) {
			if (exportToText() == true) {
				saveButton.setEnabled(false);
			}
		} else if (source == closeButton) {
			frame.dispose();
		}
	}

	public void caretUpdate(CaretEvent e) {
		saveButton.setEnabled(true);
	}

	public boolean hasFrame() {
		return frame != null;
	}

	public JFrame getFrame() {
		return frame;
	}
}
