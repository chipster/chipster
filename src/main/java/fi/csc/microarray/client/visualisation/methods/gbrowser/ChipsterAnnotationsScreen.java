package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.Time;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXTable;
import org.jdesktop.swingx.hyperlink.LinkModel;
import org.jdesktop.swingx.hyperlink.LinkModelAction;
import org.jdesktop.swingx.renderer.DefaultTableRenderer;
import org.jdesktop.swingx.renderer.HyperlinkProvider;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.screen.ScreenBase;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager;

/**
 * Screen for downloading annotations files used by the genome browser.
 * 
 * TODO Currently not used, Chipster dependencies should be removed somehow.
 * 
 * @author Taavi Hupponen
 * 
 */
public class ChipsterAnnotationsScreen extends ScreenBase implements ActionListener,
		ListSelectionListener {

	private static Logger logger = Logger.getLogger(ChipsterAnnotationsScreen.class);

	private Dimension BUTTON_SIZE = new Dimension(120, 22);
	private JFrame frame = new JFrame("Annotations");

	private JTextArea detailsTextArea = new JTextArea();
	private JScrollPane detailsScroller;

	private JButton detailsButton;
	private JButton closeButton;

	private JLabel operationLabel = new JLabel(" ");
	private JLabel parametersLabel = new JLabel(" ");
	private JLabel statusLabel = new JLabel(" ");
	private JLabel timeLabel = new JLabel(" ");
	private JLabel infoLabel = new JLabel(" ");

	private JXTable table;
	private AnnotationsTableModel tableModel;

	private TaskExecutor taskExecutor;
	private List<Task> tasks = new LinkedList<Task>();
	private AnnotationManager annotations;
	
	private enum Column {
		GENOME("Genome"), ANNOTATIONS("Annotations"), REFERENCE("Reference sequence");
		private String asString;

		Column(String asString) {
			this.asString = asString;
		}

		public String toString() {
			return asString;
		}
	};

	private class AnnotationsTableModel implements TableModel {

		private LinkedList<TableModelListener> listeners = new LinkedList<TableModelListener>();

		public boolean isCellEditable(int row, int column) {
			return false;
		}

		@Override
		public int getRowCount() {
			return annotations.getGenomes().size();
		}

		public Object getValueAt(int row, int column) {

			// Too big or negative row number
			if (row >= getRowCount() || row < 0) {
				return null;
			} 

			Column col = Column.values()[column];
			
			if (col == Column.GENOME) {
				return annotations.getGenomes().get(row).speciesId + " " + annotations.getGenomes().get(row).versionId;

			} else if (col == Column.ANNOTATIONS) {
				if (annotations.hasLocalAnnotations(annotations.getGenomes().get(row))) {
					return new LinkModel("local");
				} else {
					return new LinkModel("Download");
				}

			} else if (col == Column.REFERENCE) {
				return "local";
			} else {
				throw new IllegalArgumentException("illegal column " + column);
			}
		}

		public void addTableModelListener(TableModelListener l) {
			listeners.add(l);
		}

		public Class<?> getColumnClass(int columnIndex) {
			return getValueAt(0, columnIndex).getClass();
		}

		public int getColumnCount() {
			return Column.values().length;
		}

		public String getColumnName(int columnIndex) {
			return Column.values()[columnIndex].toString();
		}

		public void removeTableModelListener(TableModelListener l) {
			listeners.remove(l);
		}

		public void notifyListeners() {
			for (TableModelListener listener : listeners) {
				TableModelEvent tableModelEvent = new TableModelEvent(this);
				listener.tableChanged(tableModelEvent);
			}
		}

		public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
		}
	}

	@SuppressWarnings("serial")
	public ChipsterAnnotationsScreen(AnnotationManager annotations) {

		
		SwingClientApplication.setPlastic3DLookAndFeel(frame);
		frame.setPreferredSize(new Dimension(640, 480));

		frame.setLocationByPlatform(true);

		this.annotations = annotations;
		
		table = this.getTable();

		// Blue selection causes lot of visual problems, and there isn't
		// any meening for the selection
		table.setSelectionBackground(table.getBackground());
		table.setSelectionForeground(table.getForeground());

		JScrollPane tableScroller = new JScrollPane(table);
		closeButton = new JButton("Close");
		closeButton.addActionListener(this);
		closeButton.setPreferredSize(BUTTON_SIZE);

		detailsScroller = new JScrollPane(detailsTextArea);
		detailsButton = new JButton("Show Details");

		// to make sure that details get enough area
		detailsScroller.setMinimumSize(new Dimension(0, 200));
		detailsScroller.setVisible(false);
		detailsTextArea.setEditable(false);

		detailsButton.addActionListener(this);

		detailsButton.setPreferredSize(BUTTON_SIZE);

		detailsButton.setEnabled(false);

		frame.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.NORTHWEST;
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 4;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0;
		c.weighty = 1.0;
		frame.add(tableScroller, c);

		c.gridy++;
		// c.gridx=3;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.weightx = 0.0;
		c.weighty = 0.0;
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.LINE_END;
		c.insets.bottom = 8;
		c.insets.right = 8;
		c.insets.top = 8;

		frame.add(closeButton, c);

		frame.pack();
	}

	@SuppressWarnings("serial")
	private JXTable getTable() {
		if (table == null) {
			this.tableModel = new AnnotationsTableModel();
			this.table = new JXTable(tableModel) {

				// To make the custom selection colors also after font size
				// change
				@Override
				public void updateUI() {
					super.updateUI();
					setSelectionBackground(getBackground());
					setSelectionForeground(getForeground());
				}

				public Component prepareRenderer(TableCellRenderer renderer,
						int rowIndex, int vColIndex) {
					Component c = super.prepareRenderer(renderer, rowIndex,
							vColIndex);

//					Column col = Column.values()[table
//							.convertColumnIndexToModel(vColIndex)];

					rowIndex = table.convertRowIndexToModel(rowIndex);

//					if (c instanceof JComponent) {
//						JComponent jc = (JComponent) c;
//
//						Task task = tasks.get(rowIndex);
//						if (col == Column.TOOL) {
//							try {
//								jc.setToolTipText(task.getParameters()
//										.toString());
//							} catch (Exception e) {
//							}
//						} else if (col == Column.STATUS) {
//							String status = task.getState().toString();
//
//							if (task.getStateDetail() != null
//									&& task.getStateDetail().length() > 0) {
//								status += " ( " + task.getStateDetail() + " )";
//							} else if (tasks.get(rowIndex)
//									.getCompletionPercentage() != -1) {
//								status += " ( "
//										+ task.getCompletionPercentage()
//										+ "% )";
//							}
//							jc.setToolTipText(status);
//
//						} else if (col == Column.TIME) {
//							long longTime = task.getExecutionTime();
//							String min = Strings.toString(
//									(int) (longTime / 1000) / 60, 2);
//							String sec = Strings.toString(
//									(int) (longTime / 1000) % 60, 2);
//							jc.setToolTipText("Execution time: " + min + ":"
//									+ sec);
//						}
//					}
					return c;
				}
			};

			table.getSelectionModel().addListSelectionListener(this);
			table.getColumnModel().getSelectionModel()
					.addListSelectionListener(this);
			table.getSelectionModel().setSelectionMode(
					ListSelectionModel.SINGLE_SELECTION);
			table.setAutoResizeMode(JXTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
			table.getColumnModel().getColumn(Column.GENOME.ordinal())
					.setPreferredWidth(250);
			table.getColumnModel().getColumn(Column.ANNOTATIONS.ordinal())
					.setPreferredWidth(80);
			table.getColumnModel().getColumn(Column.REFERENCE.ordinal())
					.setPreferredWidth(80);
			table.setShowVerticalLines(false);

			LinkModelAction<LinkModel> linkAction = new LinkModelAction<LinkModel>() {
				public void actionPerformed(ActionEvent e) {
					logger.debug("Canceling task: "
							+ tasks.get(table.convertRowIndexToModel(table
									.getSelectedRow())));
					taskExecutor.kill(tasks.get(table
							.convertRowIndexToModel(table.getSelectedRow())));
					setVisited(true);
				}
			};

			table.getColumn(Column.ANNOTATIONS.ordinal())
					.setCellRenderer(
							new DefaultTableRenderer(new HyperlinkProvider(
									linkAction)));
		}
		return table;
	}

	private void refreshLabels(Task task) {
		detailsTextArea.setText("");
		if (task == null) {
			operationLabel.setText(" ");
			parametersLabel.setText(" ");
			statusLabel.setText(" ");
			timeLabel.setText(" ");
			infoLabel.setText(" ");
		} else {
			operationLabel.setText(task.getName());
			try {
				// TODO Task should give the List of Parameter objects instead
				// of Strings to show also the names of the parameters, like in
				// DetailsPanel
				parametersLabel.setText(task.getParameters().toString());
			} catch (Exception e) {
				parametersLabel.setText("?");
			}
			statusLabel.setText(task.getState().toString());
			timeLabel.setText((new Time(task.getStartTime())).toString());
			infoLabel.setText(task.getStateDetail());
			detailsTextArea.setText(task.getScreenOutput());

			// button is enabled if there is something to show, hiding details
			// is always possible
			detailsButton.setEnabled((task.getScreenOutput() != null && !task
					.getScreenOutput().equals(""))
					|| detailsScroller.isVisible());
		}
	}

	/**
	 * This is called when there are changes in the task information (new task,
	 * or status change) by SwingClientApplication. This should make the UI to
	 * show the changes. Running tasks that aren't already in the table are
	 * added there.
	 */
	public void refreshTasks() {
		logger.debug("Refreshing tasks in Task manager");
		for (Task task : taskExecutor.getTasks(true, false)) {
			if (!tasks.contains(task)) {
				logger.debug("\tNew task added: " + task.getName());
				tasks.add(task);
			}
		}
		tableModel.notifyListeners();
		logger.debug("Refreshing done");
	}

	public boolean hasFrame() {
		return frame != null;
	}

	public JFrame getFrame() {
		return frame;
	}

	public void actionPerformed(ActionEvent e) {

		if (e.getSource() == detailsButton) {

			if (detailsScroller.isVisible()) {
				detailsScroller.setVisible(false);
				detailsButton.setText("Show details");
			} else {
				detailsScroller.setVisible(true);
				detailsButton.setText("Hide details");
			}
		} else if (e.getSource() == closeButton) {
			frame.dispose();
		}
	}

	public void valueChanged(ListSelectionEvent e) {
		if (table.getSelectedRow() >= 0
				&& table.getSelectedRow() < tasks.size()) {
			refreshLabels(tasks.get(table.convertRowIndexToModel(table
					.getSelectedRow())));
		} else {
			refreshLabels(null);
		}
	}

	/**
	 * To be used for additional information in statusbar
	 * 
	 * @return number of failed tasks in task manager
	 */
	public int getFailedCount() {
		int i = 0;
		for (Task task : tasks) {
			if (task.getState() == Task.State.FAILED
					|| task.getState() == Task.State.TIMEOUT
					|| task.getState() == Task.State.FAILED_USER_ERROR) {
				i++;
			}
		}
		return i;
	}
}
