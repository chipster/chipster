package fi.csc.microarray.wizard.affymetrix;

import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

public class WaitDialog extends JDialog {
	private JProgressBar statusIndicator = new JProgressBar(JProgressBar.HORIZONTAL);
	private JPanel upperPanel = new JPanel();
	private JPanel lowerPanel = new JPanel();
	private JButton closeButton = new JButton("Close");
	
	public WaitDialog(JFrame frame) {
		super(frame, "Job status", true);
		
		this.getContentPane().add(lowerPanel, "South");
		this.getContentPane().add(upperPanel, "Center");
		
		upperPanel.setBorder(new EmptyBorder(new Insets(15, 10, 5, 10)));
		lowerPanel.setBorder(new EmptyBorder(new Insets(5, 10, 10, 10)));
		
		statusIndicator.setStringPainted(true);
		statusIndicator.setIndeterminate(true);
		statusIndicator.setString("Job is running");
		upperPanel.add(statusIndicator);
		
		closeButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				dispose();
			}
		});
		closeButton.setEnabled(true);
		lowerPanel.add(closeButton);
		
		setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		setModal(false);
		pack();
		setResizable(false);
		setLocationRelativeTo(frame);
		setVisible(true);
	}
	
	public void setReady() {
		statusIndicator.setIndeterminate(false);
		statusIndicator.setString("Ready");
		closeButton.setEnabled(true);
	}
}
