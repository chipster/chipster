package fi.csc.microarray.client.dialog;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dialog.DialogInfo.Type;

public class ChipsterDialog extends JDialog {
	
	public enum DetailsVisibility { 
		DETAILS_ALWAYS_VISIBLE,
		DETAILS_ALWAYS_HIDDEN,
		DETAILS_VISIBLE,
		DETAILS_HIDDEN;
		
		public boolean isButtonEnabled() {
			return this == DETAILS_HIDDEN || this == DETAILS_VISIBLE;  
		}
		
		public boolean isInitiallyVisible() {
			return this == DETAILS_ALWAYS_VISIBLE || this == DETAILS_VISIBLE;
		}
	}; 
	
	
	private static final int DETAILS_AREA_HEIGHT = 200;
	private static final int DETAIL_AREA_WIDTH = 400;

	private JPanel detailsPanel = new JPanel();
	private SwingClientApplication application;
	
	private DetailsVisibility detailsVisibility;
	private DialogCloseListener dialogCloseListener;
	
	public ChipsterDialog(SwingClientApplication app, DialogInfo dialogInfo, 
	        DetailsVisibility detailsVisibility) {
	    
        super(app != null ? app.getMainFrame() : null);
	    this.application = app;
		this.detailsVisibility = detailsVisibility;
		
		// initialise dialog layout
		SwingClientApplication.setPlastic3DLookAndFeel(this);
		
		// insert components to layout
		JPanel mainPanel = getMainPanel(dialogInfo);
		mainPanel.setBorder(new EmptyBorder(10, 10, 0, 0));
		this.getContentPane().add(mainPanel);
		
		// locate and pack
		this.setLocationByPlatform(true);
		this.pack();
	}

	private JPanel getMainPanel(DialogInfo dialogInfo) {
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new GridBagLayout());
		GridBagConstraints g = new GridBagConstraints();
		g.anchor = GridBagConstraints.NORTHWEST;
		g.gridx = 0;
		g.gridy = 0;
		g.weightx = 0.0;
		g.insets = new Insets(5, 0, 10, 5);
		Icon icon;
		switch (dialogInfo.getSeverity()) {
		case INFO:
			icon = UIManager.getLookAndFeelDefaults().getIcon("OptionPane.informationIcon");
			break;
		case WARNING:
			icon = UIManager.getLookAndFeelDefaults().getIcon("OptionPane.warningIcon");
			break;
		case ERROR:
			icon = UIManager.getLookAndFeelDefaults().getIcon("OptionPane.errorIcon");
			break;
		case QUESTION:
			icon = UIManager.getLookAndFeelDefaults().getIcon("OptionPane.questionIcon");
		default:
			throw new IllegalArgumentException("unsupported severity level: " + dialogInfo.getSeverity());
		}
		mainPanel.add(new JLabel(icon), g);
		g.gridx++;
		g.weightx = 1.0;
		g.gridwidth = 2;
		JTextArea titleArea = new JTextArea(dialogInfo.getTitle());
		titleArea.setEditable(false);
		titleArea.setLineWrap(true);
		titleArea.setWrapStyleWord(true);
		titleArea.setFont(titleArea.getFont().deriveFont(Font.BOLD, 13f));
		titleArea.setColumns(32);
		titleArea.setBackground(mainPanel.getBackground());

		JTextArea messageArea = new JTextArea(dialogInfo.getMessage());
		messageArea.setEditable(false);
		messageArea.setLineWrap(true);
		messageArea.setWrapStyleWord(true);
		messageArea.setColumns(38);
		messageArea.setBackground(mainPanel.getBackground());
		
		mainPanel.add(titleArea, g);
		g.gridy++;
		mainPanel.add(messageArea, g);
		g.gridy++;
		g.gridwidth = 1;
		g.gridx = 1;
		g.weightx = 0.0;
		g.insets = new Insets(5, 0, 0, 5);
		
		if (detailsVisibility == DetailsVisibility.DETAILS_VISIBLE ||
				detailsVisibility == DetailsVisibility.DETAILS_HIDDEN) {
			mainPanel.add(getDetailsButton(), g);
		}
		
        g.gridx++;
		if (dialogInfo.getFeedbackVisible()) {
		    JButton feedbackButton = new JButton("Send a report");
		    feedbackButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    FeedbackDialog feedback = new FeedbackDialog(application);
                    feedback.showDialog();
                }
            });
		    mainPanel.add(feedbackButton, g);
		}
		
        g.gridx++;
		if (dialogInfo.getType() == Type.OPTION) {
			g.anchor = GridBagConstraints.EAST;
			JButton cancelButton = new JButton("Cancel");
			cancelButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ChipsterDialog.this.dispose();
					if (dialogCloseListener != null) {
						dialogCloseListener.dialogClosed(false);
					}
				}
			});
			mainPanel.add(cancelButton, g);
			g.anchor = GridBagConstraints.WEST;
		}

		g.gridx++;
		String okButtonContent = dialogInfo.getType().getButtonText();
		if (okButtonContent != null) {
			JButton okButton = new JButton(dialogInfo.getType().getButtonText());
			okButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ChipsterDialog.this.dispose();
					if (dialogCloseListener != null) {
						dialogCloseListener.dialogClosed(true);
					}
				}
			});
			mainPanel.add(okButton, g);
		}
		g.gridy++;
		g.gridx = 1;
		g.fill = GridBagConstraints.BOTH;
		g.gridwidth = 4;
		g.weighty = 1.0;
		
		detailsPanel.setLayout(new CardLayout());
		
		JTextArea detailsArea = new JTextArea(dialogInfo.getDetails());		
		detailsArea.setBorder(new LineBorder(Color.BLACK));
		detailsArea.setEditable(false);
		JScrollPane scrollPane = new JScrollPane(detailsArea);
		detailsPanel.add(scrollPane, "visible");
		
		JPanel placeholder = new JPanel();
		detailsPanel.add(placeholder, "hidden");
		
		setDetailsVisible(false);
		mainPanel.add(detailsPanel, g);

		return mainPanel;
	}

	private JButton getDetailsButton() {
		final JButton button;
		if (detailsVisibility == DetailsVisibility.DETAILS_HIDDEN) {
			button = new JButton("Show details");
		} else {
			button = new JButton("Hide details");
		}
		
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (button.getText().equals("Show details")) {
					setDetailsVisible(true);
					button.setText("Hide details");
				} else {
					setDetailsVisible(false);
					button.setText("Show details");					
				}
			}
		});
		return button;
	}
	
	public void setDetailsVisible(boolean visible) {
	    CardLayout cl = (CardLayout)(detailsPanel.getLayout());
	    cl.show(detailsPanel, visible ? "visible" : "hidden");
	    detailsPanel.setPreferredSize(new Dimension(DETAIL_AREA_WIDTH, visible ? DETAILS_AREA_HEIGHT : 0));
	    this.pack();
	}
	
	public static void showDialog(SwingClientApplication application,
	        DialogInfo dialogInfo, DetailsVisibility detailsVisibility,
	        boolean modal) {
		showDialog(application, dialogInfo, detailsVisibility, modal, null);
	}
	
	public static void showDialog(SwingClientApplication application,
	        DialogInfo dialogInfo, DetailsVisibility detailsVisibility,
	        boolean modal, DialogCloseListener dialogCloseListener) {
	    
		ChipsterDialog dialog = new ChipsterDialog(application, dialogInfo, detailsVisibility);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.setModal(modal);
		dialog.setDialogCloseListener(dialogCloseListener);
		
		dialog.setDetailsVisible(detailsVisibility.isInitiallyVisible());
		dialog.getDetailsButton().setEnabled(detailsVisibility.isButtonEnabled());
		
		dialog.setVisible(true);	
	}

	public static interface DialogCloseListener {
		public void dialogClosed(boolean okSelected);
	}
	
	public void setDialogCloseListener(DialogCloseListener dialogCloseListener) {
		this.dialogCloseListener = dialogCloseListener;
	}	
}
