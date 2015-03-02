package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.Timer;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.SystemMonitorUtil;

public class StatusBar {

	private SwingClientApplication application;
	private JPanel statusPanel = null;
	private JLabel statusLabel = null;
	private JProgressBar jobStatusIndicator = null;
	private Color normalStatusIndicatorColor;
	private JProgressBar memoryIndicator;
	private JButton jobListButton;

	private int oldTaskCount = 0;
	private Timer blinker;

	
	public StatusBar(SwingClientApplication application) {
		this.application = application;
	}
	
	public JPanel getStatusPanel() {
		if (statusPanel == null) {
			// add status bar
			jobStatusIndicator = new JProgressBar();

			jobStatusIndicator.setName("jobStatusIndicator");
			jobStatusIndicator.setStringPainted(true);
			normalStatusIndicatorColor = jobStatusIndicator.getForeground();
			
			memoryIndicator = new JProgressBar(0, 100);
			memoryIndicator.setStringPainted(true);
			memoryIndicator.setToolTipText("Shows how much of the available memory is allocated. Click the indicator to release unused memory (garbage collection).");
			memoryIndicator.addMouseListener(new MouseListener() {
				public void mouseClicked(MouseEvent e) {
					application.garbageCollect();
				}

				public void mouseEntered(MouseEvent e) { /* do nothing */
				}

				public void mouseExited(MouseEvent e) { /* do nothing */
				}

				public void mousePressed(MouseEvent e) { /* do nothing */
				}

				public void mouseReleased(MouseEvent e) { /* do nothing */
				}
			});

			statusLabel = new JLabel();
			String labelText = "Connected to " + DirectoryLayout.getInstance().getConfiguration().getString("messaging", "broker-host");
			statusLabel.setText(labelText);
			statusLabel.setBorder(jobStatusIndicator.getBorder());

			jobListButton = new JButton("View jobs");
			jobListButton.setName("jobListButton");
			jobListButton.setToolTipText("View jobs");

			jobListButton.addMouseListener(new MouseListener() {

				public void mouseClicked(MouseEvent e) {
					viewTasks();
				}

				public void mouseEntered(MouseEvent e) {
				}

				public void mouseExited(MouseEvent e) {
				}

				public void mousePressed(MouseEvent e) {
				}

				public void mouseReleased(MouseEvent e) {
				}
			});

			this.taskCountChanged(0, 0, false);

			statusPanel = new JPanel(new GridBagLayout());
			GridBagConstraints c = new GridBagConstraints();
			c.weightx = 1.0;
			c.weighty = 1.0;
			c.fill = GridBagConstraints.BOTH;
			c.insets.set(0, 0, 0, 5);
			statusPanel.add(statusLabel, c);
			c.gridx = 1;
			c.weightx = 0.0;
			statusPanel.add(jobListButton, c);
			c.gridx = 2;
			statusPanel.add(jobStatusIndicator, c);
			c.insets.set(0, 0, 0, 0);
			c.gridx = 3;
			memoryIndicator.setPreferredSize(new Dimension(200, 20));
			statusPanel.add(memoryIndicator, c);
		}
		return statusPanel;
	}

	public void viewTasks() {
		application.refreshTaskList();

		// show
		application.getTaskListScreen().getFrame().setVisible(true);
		application.getTaskListScreen().getFrame().setExtendedState(JFrame.NORMAL);
		application.getTaskListScreen().getFrame().toFront();
		application.getTaskListScreen().getFrame().setFocusable(true);
		application.getTaskListScreen().getFrame().requestFocus();
		application.getTaskListScreen().getFrame().toFront();
		jobListButton.setToolTipText("Hide Task manager");
	}
	
	private void setProgressBar(String indicatorText, boolean indeterminate, int value) {
		jobStatusIndicator.setString(indicatorText);
		jobStatusIndicator.setIndeterminate(indeterminate);
		jobStatusIndicator.setValue(value);
		if (value > 0) {
			jobStatusIndicator.setForeground(normalStatusIndicatorColor.darker());
		} else {
			jobStatusIndicator.setForeground(normalStatusIndicatorColor);
		}
	}
	
	protected void taskCountChanged(int taskCount, int completion, boolean useBlinking) {

		// update status bar
		boolean jobStarted;
		if (oldTaskCount < taskCount) {
			jobStarted = true;
		} else {
			jobStarted = false;
		}

		oldTaskCount = taskCount;

		// update progress bar state
		if (taskCount == 0) {
			setProgressBar("0 jobs running", false, 0);
		} else {
			if (taskCount == 1) {
				if (completion != -1) {
					setProgressBar("Transferring data", false, completion);
				} else {
					setProgressBar(taskCount + " job running", true, 0);
				}
			} else {
				setProgressBar(taskCount + " jobs running", true, 0);
			}
			
		}		

		if (useBlinking && jobStarted) {
			blinkProgressBar();
		}

		// update title
		application.updateWindowTitleJobCount(taskCount);

		// update task list
		application.refreshTaskList();
	}

	/**
	 * Adds a fancy blinking effect to progress bar after job is started
	 * 
	 */
	private void blinkProgressBar() {

		class Blinker implements ActionListener {

			private boolean isBold;
			private int repeats;
			private int repeated;

			public Blinker(int repeats) {
				this.isBold = false;
				this.repeats = repeats;
				this.repeated = 0;
			}

			public void actionPerformed(ActionEvent e) {
				if (repeated == repeats) {
					jobStatusIndicator.setFont(jobStatusIndicator.getFont().deriveFont(Font.PLAIN));
					((Timer) e.getSource()).stop();
					return;
				}

				if (isBold) {
					jobStatusIndicator.setFont(jobStatusIndicator.getFont().deriveFont(Font.PLAIN));
					isBold = false;
				} else {
					jobStatusIndicator.setFont(jobStatusIndicator.getFont().deriveFont(Font.BOLD));
					isBold = true;
				}
				repeated++;
			}
		}

		if (blinker != null) {
			blinker.stop();
			jobStatusIndicator.setFont(jobStatusIndicator.getFont().deriveFont(Font.PLAIN));
		}

		blinker = new Timer(350, new Blinker(6));
		blinker.start();

	}

	public void updateMemoryIndicator() {
		memoryIndicator.setString("Used memory " + SystemMonitorUtil.getMemInfo());
		memoryIndicator.setValue((int) (((float) SystemMonitorUtil.getUsed()) / ((float) Runtime.getRuntime().maxMemory()) * 100f));
	}

	public void setFontSize(float fontSize) {
		jobStatusIndicator.setFont(jobStatusIndicator.getFont().deriveFont(fontSize));
	}

}
