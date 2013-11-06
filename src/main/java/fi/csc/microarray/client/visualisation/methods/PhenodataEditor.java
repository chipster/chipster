package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.CellEditorListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.visualisation.ExtendedJXTable;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.features.table.TableBeanEditor;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule;


/**
 * Editor for phenodata type of tabular DataBeans.
 * 
 * @author Aleksi Kallio, Petri Klemelä
 *
 */
public class PhenodataEditor extends Visualisation implements DataChangeListener {
	
	public static final String PHENODATA_SAMPLE_COLUMN = "sample";
	public static final String PHENODATA_NAME_COLUMN = "original_name";
	public static final String PHENODATA_DESCRIPTION_COLUMN = "description";
	private static final String PHENODATA_CHIPTYPE_COLUMN = "chiptype";
	private static final String PHENODATA_GROUP_COLUMN = "group";
	
	public static boolean isEditablePhenodataColumn(String columnName) {
		return !PHENODATA_SAMPLE_COLUMN.equals(columnName) && !PHENODATA_NAME_COLUMN.equals(columnName);
	}
	
	public static boolean isGroupPhenodataColumn(String columnName) {
		return PHENODATA_GROUP_COLUMN.equals(columnName);
	}

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	private JPanel paramPanel;
	private JTextField columnNameField;
	private JComboBox columnsCombo;
	
	
	
	@SuppressWarnings("serial")
    private class PhenodataTable extends ExtendedJXTable {
		private static final int NO_SCROLL_WIDTH = 400;
		private static final int IDENTIFIER_COLUMN_WIDTH = 100;
				
		private int lastX = -1;
		private int lastY = -1;

		public PhenodataTable(DataBean data) {
			super(data);			
			this.setHorizontalScrollEnabled(this.getWidth() < NO_SCROLL_WIDTH);			
			
			this.addKeyListener(new KeyListener(){
				public void keyPressed(KeyEvent e) {
				}

				public void keyReleased(KeyEvent e) {
				}

				public void keyTyped(KeyEvent e) {										

					if(getCellEditor() != null){						
						
						int y = getEditingRow();
						int x = getEditingColumn();
						
						//Cell in sorted column can't be emptied, because this would 
						//move the emptied cell to different place.
						if(getSortedColumn() != null && getSortedColumn().equals(getColumn(x))){
							return;
						}
						
						if(y != lastY || x != lastX){
														
							getCellEditor().stopCellEditing();
												
							lastY = y;
							lastX = x;
														
							setValueAt("", y, x);
							
							boolean success = editCellAt(y, x);
							
							if (success) {
								changeSelection(y, x, false, false);
							}
						}
						
						//May be null after stopCellEditing
						if(getCellEditor() != null){
							getCellEditor().addCellEditorListener(new CellEditorListener(){

								public void editingCanceled(ChangeEvent e) {
									lastX = lastY = -1;
								}

								public void editingStopped(ChangeEvent e) {
									lastX = lastY = -1;
								}
							});
						}
					}
				}				
			});
			
		}
		
		@Override
		public void setModel(TableModel model){
			super.setModel(model);
			setColumnsComboContent();
		}

		private Color uneditableBg = getBackground();
		
		public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
		    
			// Returns cell under the event location
			Component cell = super.prepareRenderer(renderer, row, column);
			
		    if (this.isColumnSelected(column) && this.isRowSelected(row)){
		        // Cell is selected
		    	cell.setBackground(this.getSelectionBackground());
		    	cell.setForeground(this.getSelectionForeground());
		    } else {
		    	if(!isEditablePhenodataColumn(getColumnName(column))){
		    		cell.setBackground(this.getUneditableBackground());
		    	} else {
		    		cell.setBackground(this.getBackground());
		    	}
		        cell.setForeground(this.getForeground());
		    }
		    return cell;
		}		
		
		public Color getUneditableBackground(){
			return this.uneditableBg;
		}
		
		public void setUneditableBackground(Color c){
			this.uneditableBg = c;
		}
	}
	
	private void setColumnsComboContent(){
		
		columnsCombo.removeAllItems();
		Iterator<String> iter = tableEditor.getEditable().getColumnNames().iterator();
		while (iter.hasNext()) {
			String columnName = iter.next();
			if (isEditablePhenodataColumn(columnName)) {
				columnsCombo.addItem(columnName);
			}
		}
	}
	
	@Override
	public JPanel getParameterPanel() {
		logger.debug("in PhenodataEditor.getParameterPanel()");
		if (paramPanel == null) {
			logger.debug("in PhenodataEditor.getParameterPanel() paramPanel is null");
			paramPanel =  new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			
			columnNameField = new JTextField("new_column");
			columnsCombo = new JComboBox();
			
			JButton addButton = new JButton("Add");
			addButton.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					
					if(columnNameField.getText() == null || 
							columnNameField.getText().equals("") ||							
							tableEditor.getEditable().containsColumn(columnNameField.getText())) {
						JOptionPane.showMessageDialog(
								((SwingClientApplication)application).getMainFrame(), 
								"Name of the column is empty or exists already.");
						return;
					}
					try {
						tableModel.addColumn(columnNameField.getText());
					} catch (Exception e1) {
						application.reportException(e1);
					}
					
					table.createDefaultColumnsFromModel();
					updatePhenodataTableHeaders();					
					setColumnsComboContent();
				}
			});
			
			JButton removeButton = new JButton("Remove");
			removeButton.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {

					try {
						tableModel.removeColumn((String)columnsCombo.getSelectedItem());
					} catch (Exception e1) {
						application.reportException(e1);
					}
					
					table.createDefaultColumnsFromModel();
					updatePhenodataTableHeaders();
					setColumnsComboContent();
				}
			});

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.insets.set(5, 10, 5, 10);
			c.anchor = GridBagConstraints.NORTHWEST;
			c.weightx = 1.0;
			c.fill = GridBagConstraints.HORIZONTAL;
			paramPanel.add(new JLabel("Add a new column: "),c);
			c.gridy++;
			paramPanel.add(columnNameField,c);
			c.gridy++;
			paramPanel.add(addButton,c);
			c.gridy++;
			paramPanel.add(new JLabel("Remove column: "),c);
			c.gridy++;
			paramPanel.add(columnsCombo,c);
			c.gridy++;
			paramPanel.add(removeButton,c);
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;				
			paramPanel.add(new JPanel(),c);
			logger.debug("in PhenodataEditor.getParameterPanel() components created and added");
		}
		return paramPanel;
	}
	
	/**
	 * PopupMenu for phenodata editor panel. Allows copy and paste operations.
	 * 
	 * @author Mikko Koski
	 *
	 */
	@SuppressWarnings("serial")
    public class PhenodataPopupMenu extends JPopupMenu implements ActionListener {

		private ExtendedJXTable table;
		private JMenuItem copyMenuItem;
		private JMenuItem pasteMenuItem;
		
		public PhenodataPopupMenu(ExtendedJXTable table) {
			this.table = table;
			copyMenuItem = new JMenuItem("Copy");
			pasteMenuItem = new JMenuItem("Paste");
			copyMenuItem.addActionListener(this);
			pasteMenuItem.addActionListener(this);
			
			this.add(copyMenuItem);
			this.add(pasteMenuItem);
		}
		
		public void actionPerformed(ActionEvent e) {
			if (e.getSource() == copyMenuItem) {
				table.copy();
			}
			if (e.getSource() == pasteMenuItem) {
				table.paste();
			}
		}
	}
	
	private final class PhenodataTableModel implements TableModel {

		private final TableBeanEditor tableEditor;
		private LinkedList<TableModelListener> listeners = new LinkedList<TableModelListener>();
		
		private PhenodataTableModel(TableBeanEditor tableEditor) {
			super();
			this.tableEditor = tableEditor;				
		}

		public int getRowCount() {
			return tableEditor.getEditable().getRowCount();
		}

		public int getColumnCount() {
			return tableEditor.getEditable().getColumnCount();
		}

		public Class<?> getColumnClass(int columnIndex) {
			return Object.class;
		}

		public boolean isCellEditable(int rowIndex, int columnIndex) {
			return isEditablePhenodataColumn(tableEditor.getEditable().getColumnName(columnIndex));
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			String title = tableEditor.getEditable().getColumnName(columnIndex);
			return tableEditor.getEditable().getValue(title, rowIndex);
		}

		public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
			try {
				String title = tableEditor.getEditable().getColumnName(columnIndex);
				logger.debug("set column " + title + " at row " + rowIndex + " to value " + aValue);
				tableEditor.getEditable().setValue(title, rowIndex, (String)aValue);
				tableEditor.write();
				for (TableModelListener listener : listeners) {
					listener.tableChanged(new TableModelEvent(this));
				}
			} catch (Exception e) {
				application.reportException(e);
			}
		}

		public void addTableModelListener(TableModelListener l) {
			listeners.add(l);
		}

		public void removeTableModelListener(TableModelListener l) {
			listeners.remove(l);
		}
		
		public void addColumn(String title) throws Exception{
			List<String> values = new ArrayList<String>();
			for ( int i = 0; i < tableEditor.getEditable().getRowCount() ; i++){
				values.add("");
			}
			tableEditor.getEditable().addColumn(title, values);
			tableEditor.write();
		}
		
		public void removeColumn(String title) throws Exception{
			tableEditor.getEditable().removeColumn((String)columnsCombo.getSelectedItem());			
			tableEditor.write();
		}

		public String getColumnName(int columnIndex) {
			return tableEditor.getEditable().getColumnName(columnIndex);
		}
	}

	private static final Logger logger = Logger.getLogger(PhenodataEditor.class);

	private DataBean data;
	private PhenodataTable table;
	private TableBeanEditor tableEditor;
	private PhenodataTableModel tableModel;
	
	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
			
		// initialise data
		this.data = data;
		this.tableEditor = new TableBeanEditor(data);
		this.tableModel = new PhenodataTableModel(tableEditor);
		
		// initialise UI
		this.table = new PhenodataTable(data);
		//Dimension tableSize = size;
		//tableSize.setSize(tableSize.getWidth()-MARGIN, tableSize.getHeight()-MARGIN);
		table.setUneditableBackground(VisualConstants.PHENODATA_TABLE_UNEDITABLE_CELL_BACKGROUND);
	
		table.setModel(tableModel);
		table.getColumn(0).setPreferredWidth(PhenodataTable.IDENTIFIER_COLUMN_WIDTH);		
		
		// set warning signs
		updatePhenodataTableHeaders();
		
		// listen for phenodata content change
		application.getDataManager().addDataChangeListener(this);
		
		//table.setPreferredScrollableViewportSize(tableSize);
		table.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		table.setColumnControlVisible(true);
		JScrollPane tableScroller = new JScrollPane(table);
        
        // adds popupmenu listener
		table.addMouseListener(new MouseAdapter(){
			
			@Override
	        public void mousePressed(MouseEvent e) {
	        	maybeShowPopup(e);        
	        }
	        
			@Override
	        public void mouseReleased(MouseEvent e) {
	        	maybeShowPopup(e);        
	        }
	            		    				
		    private void maybeShowPopup(MouseEvent e) {
		        if (e.isPopupTrigger()) {
		        	JPopupMenu popup = new PhenodataPopupMenu(table);
	                popup.show(e.getComponent(), e.getX(), e.getY());
		        }
		    }
		});										
		
		return tableScroller;
	}
	
	private void updatePhenodataTableHeaders(){
		// Set warning icon
		for (Object columnObject : table.getColumns()){
			if (columnObject instanceof TableColumn) {
				TableColumn tableColumn = (TableColumn) columnObject;
				if(isGroupPhenodataColumn(tableColumn.getHeaderValue().toString())){
					DefaultTableCellRenderer header = new DefaultTableCellRenderer();
					if(!data.queryFeatures("/phenodata/is-complete").exists()){
						header.setIcon(VisualConstants.PHENODATA_ICON);
						logger.debug("Header updated. Warning icon enabled.");
					} else {
						header.setIcon(null);
						logger.debug("Header updated. Warning icon disabled.");
						
					}
					
					// FIXME hackhack
					header.setBackground(VisualConstants.TEXTAREA_UNEDITABLE_BACKGROUND);
	
					header.repaint();
					table.getTableHeader().repaint();

					tableColumn.setHeaderRenderer(header);
				}
			}
		}
	}

	public void dataChanged(DataChangeEvent evt) {
		logger.debug("received " + evt.getClass().getSimpleName());
		if (evt instanceof ContentChangedEvent) {
			ContentChangedEvent cce = (ContentChangedEvent)evt;
			 if (this.data == cce.getDataItem()) {
				 updatePhenodataTableHeaders();
				 logger.debug("phenodata headers updated for " + cce.getDataItem().getName());
			 }
		}
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.hasTypeTag(BasicModule.TypeTags.PHENODATA);
	}
	

	@Override
	public void removeVisualisation(){
		application.removeClientEventListener(table);
	}
}
