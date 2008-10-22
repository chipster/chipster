package fi.csc.microarray.client.dataimport.table;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.ProgressInformator;
import fi.csc.microarray.client.dataimport.RunnableImportProcess;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

/**
 * A small dialog which implements ProcessInformator to show detailed informator 
 * of the state of the process which is running.
 * 
 * @author mkoski
 *
 */
public class InformationDialog extends JFrame implements ActionListener, ProgressInformator {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(InformationDialog.class);

	private final Dimension BUTTON_SIZE = new Dimension(80, 25);
	private final Dimension PROGRESSBAR_SIZE = new Dimension(300, 18);
	
	private RunnableImportProcess process;
	private JPanel contentPane;
	private JLabel label;
	private JButton cancelButton;
	private JProgressBar progressBar;
	
	private ImportScreen owner;
	
	public InformationDialog(String title, String message, ImportScreen owner) {
		super(title);
		this.label = new JLabel(message);
		
		// This keeps the progressbar size correct after setMessage
//		label.setPreferredSize(new Dimension(
//				(int)PROGRESSBAR_SIZE.getWidth() + 100, 
//				(int)PROGRESSBAR_SIZE.getHeight()));
		label.setPreferredSize(PROGRESSBAR_SIZE);
				
		this.owner = owner;
		
		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);
		
		progressBar = new JProgressBar();
		progressBar.setPreferredSize(PROGRESSBAR_SIZE);
		
		contentPane = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		
		// Label
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(8, 12, 10, 4);
		c.gridx = 0; 
		c.gridy = 0;
		contentPane.add(label, c);
		
		// Combobox
		c.insets.set(0, 12, 8, 12);
		c.gridy = 1;
		contentPane.add(progressBar, c);
		
		// Buttons
		c.insets.set(4, 2, 12, 8);
		c.anchor = GridBagConstraints.EAST;
		c.gridy = 2;
		contentPane.add(cancelButton, c);
		
		logger.debug("Components on the dialog: " + contentPane.getComponentCount());
		
		this.add(contentPane);
		this.pack();
		this.setLocationRelativeTo(null);
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == cancelButton){
			this.stopProcess();
		}
	}
	
	public void setValue(int value){
		class ValueChanger implements Runnable{
			private int value;
			
			public ValueChanger(int value) {
				this.value = value;
			}
			public void run() {
				progressBar.setValue(value);	
			}
		}
		SwingUtilities.invokeLater(new ValueChanger(value));
	}

	/**
	 * Sets message to the information dialog. This is a thread-safe 
	 * method
	 */
	public void setMessage(String message) {
		class MessageChanger implements Runnable{
			
			private String message;
			
			public MessageChanger(String message) {
				this.message = message;
			}

			public void run() {
				label.setText(message);
			}
		}
		SwingUtilities.invokeLater(new MessageChanger(message));
	}

	public void setMaximumValue(int max) {
		class MaximumChanger implements Runnable{
			
			private int max;
			
			public MaximumChanger(int max) {
				this.max = max;
			}
			
			public void run(){
				progressBar.setMaximum(max);
			}
		}
		SwingUtilities.invokeLater(new MaximumChanger(max));
		
	}

	public void setMinimunValue(int min) {
		class MinimumChanger implements Runnable{
			
			private int min;
			
			public MinimumChanger(int min) {
				this.min = min;
			}
			
			public void run(){
				progressBar.setMaximum(min);
			}
		}
		SwingUtilities.invokeLater(new MinimumChanger(min));
	}

	public void destroyInformator() {
		if(owner != null){
			owner.getFrame().setEnabled(true);
		}
		this.setVisible(false);
		try{
			this.dispose();
		} catch (Exception e){
			// Do nothing
		}
	}

	public void initializeInformator() {
		if(owner != null){
			owner.getFrame().setEnabled(false);
		}
		setVisible(true);

	}

	public void stopProcess() {
		if(process != null){
			this.process.stopProcess();
		} else {
			throw new IllegalStateException("No process set to this informator");
		}
	}

	public void setProcess(RunnableImportProcess process) {
		this.process = process;
	}

	public void setIndeterminate(boolean newValue) {
		this.progressBar.setIndeterminate(newValue);
	}
}
