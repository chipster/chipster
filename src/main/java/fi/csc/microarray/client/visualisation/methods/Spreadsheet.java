package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.table.DefaultTableModel;

import fi.csc.microarray.client.visualisation.MicroarrayTable;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.features.RestrictModifier;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * A GUI component for showing (and some day maybe even editing!) the
 * microarray data - that is, the channels and their intensities.
 * It basically takes the data from a dataset's microarray and puts it
 * to a JTable. Once gene names are implemented in microarrays, this
 * could be used to annotate them (write some describing names and notes
 * for different genes).
 * 
 * @author Janne KÃ¤ki, Mikko Koski, Aleksi Kallio
 *
 */


public class Spreadsheet extends Visualisation {
	
	public Spreadsheet(VisualisationFrame frame) {
		super(frame);
	}	

	/**
	 * PopupMenu for spread sheet view. Allows copy operation and annotating 
	 * using Bioconductor
	 *
	 */
	public class SpreadsheetPopupMenu extends JPopupMenu implements ActionListener {

		private MicroarrayTable table;
		private JMenuItem annotateMenuItem;
		private JMenuItem copyMenuItem;
		private JMenuItem filterMenuItem;
		
		public SpreadsheetPopupMenu(MicroarrayTable table) {
			this.table = table;
			
			copyMenuItem = new JMenuItem("Copy");			
			annotateMenuItem = new JMenuItem("Create dataset and annotate with Bioconductor");
			filterMenuItem = new JMenuItem("Create dataset");
			
			copyMenuItem.addActionListener(this);			
			annotateMenuItem.addActionListener(this);
			filterMenuItem.addActionListener(this);
			
			this.add(copyMenuItem);
			this.add(new JSeparator());
			this.add(annotateMenuItem);
			this.add(filterMenuItem);
		}
		
		public void actionPerformed(ActionEvent e) {
						
			if (e.getSource() == annotateMenuItem) {
				VisualisationUtilities.annotateBySelection(getFrame().getDatas());
			}
			if (e.getSource() == copyMenuItem) {
				table.copy();
			}
			if(e.getSource() == filterMenuItem) {
				VisualisationUtilities.filterBySelection(getFrame().getDatas());
			}
		}
	}

	private final int COLUMNS_REQUIRES_SCROLLING = 8;

	private MicroarrayTable table;
	
	/**
	 * Creates a new TablePanel, which (for now) is dataset specific.
	 * 
	 * @param datas The dataset for which to create the table.
	 * @throws Exception 
	 * @throws MicroarrayException 
	 */
	@Override
	public JComponent getVisualisation(Dataset data) throws Exception {
		JPanel panel = new JPanel(new BorderLayout());				
	
		Table rowCounter = data.queryFeatures("/column/*").asTable();
		int rowCount = 0;
		while (rowCounter.nextRow()) {
			rowCount++;
		}
		
		Table columns = data.queryFeatures("restrict(/column/*)").asTable();
		String[] columnTitles = new String[columns.getColumnCount()];
		int counter = 0;
		for (String column : columns.getColumnNames()) {
			columnTitles[counter] = column;
			counter++;
		}
		
		Object[][] rowData =
			new Object[RestrictModifier.RESTRICT_TO_ROWS < rowCount ? RestrictModifier.RESTRICT_TO_ROWS : rowCount][columns.getColumnCount()];
		int row = 0;
		while (columns.nextRow()) {
			int column = 0;
			for (String columnName : columns.getColumnNames()) {
				Object value = columns.getValue(columnName);
				if (value instanceof Float) {
					value = new NicelyShowingFloat((Float)value);
				}
				rowData[row][column] = value;
				column++;
			}
			row++;
		}
		
		table = new MicroarrayTable(data);
		DefaultTableModel tableModel = new DefaultTableModel(rowData, columnTitles) {			
			@Override
			public boolean isCellEditable(int row, int column){
				return false;
			}			
		};
		table.setModel(tableModel);
		table.setColumnControlVisible(true);
		JScrollPane tableScroller = new JScrollPane(table);
        table.setBackground(java.awt.Color.white);
		table.setHorizontalScrollEnabled(columns.getColumnCount() > COLUMNS_REQUIRES_SCROLLING);
		
		// Adds popupmenu listener
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
		        	JPopupMenu popup = new SpreadsheetPopupMenu(table);
	                popup.show(e.getComponent(), e.getX(), e.getY());
		        }
		    }
		});
		
		panel.add(new JLabel("Showing " + rowData.length + " rows of " + rowCount), BorderLayout.NORTH);
		panel.add(tableScroller, BorderLayout.CENTER);
		table.updateSelectionsFromApplication();
	
		table.sendEvents(true);
		return panel; 
	}

	
	public static class NicelyShowingFloat implements Comparable {

		private Float floatValue;
		
		public NicelyShowingFloat(Float floatValue) {
			this.floatValue = floatValue;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof NicelyShowingFloat) {
				return floatValue.equals(((NicelyShowingFloat)obj).floatValue);
			}
			throw new IllegalArgumentException("cannot compare to instance of " + obj.getClass().getSimpleName());
		}

		@Override
		public int hashCode() {
			return floatValue.hashCode();
			
		}
		
		@Override
		public String toString() {
			String s = floatValue.toString();
			
			if (s.endsWith(".0")) {
				s = s.substring(0, s.length()-2);
			}
			
			return s;
		}
		
		public int compareTo(Object obj) {
			if (obj instanceof NicelyShowingFloat) {
				return floatValue.compareTo(((NicelyShowingFloat)obj).floatValue);
			}
			throw new IllegalArgumentException("cannot compare to instance of " + obj.getClass().getSimpleName());
		}
		
	}
	
	@Override
	public boolean canVisualise(Dataset bean) throws MicroarrayException {

		if (bean.isContentTypeCompatitible("text/tab", "text/csv")) {
			return true; // clearly tabular
			
		} else if (bean.isContentTypeCompatitible("application/cel")) {
			return !bean.queryFeatures("/embedded-binary-content/").exists(); // might have embedded binary content
			
		} else {
			return false;
		}
	}
	
	@Override
	public void removeVisualisation(){
		application.removePropertyChangeListener(table);
	}
}
