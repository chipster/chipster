 package fi.csc.microarray.client.dataimport.table;

import java.awt.Color;
import java.awt.Component;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.lang.reflect.InvocationTargetException;

import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumn;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXTable;
import org.jdesktop.swingx.table.DefaultTableColumnModelExt;
import org.jdesktop.swingx.table.TableColumnModelExt;

import fi.csc.microarray.client.dataimport.ColumnType;
import fi.csc.microarray.client.dataimport.ColumnTypeManager;
import fi.csc.microarray.client.dataimport.ConversionModel;
import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.events.ColumnTitlesChangedEvent;
import fi.csc.microarray.client.dataimport.events.ConversionModelChangeListener;
import fi.csc.microarray.client.dataimport.events.DecimalSeparatorChangedEvent;
import fi.csc.microarray.client.dataimport.events.DelimiterChangedEvent;
import fi.csc.microarray.client.dataimport.events.FooterChangedEvent;
import fi.csc.microarray.client.dataimport.events.HeaderChangedEvent;
import fi.csc.microarray.client.dataimport.events.InputFileChangedEvent;
import fi.csc.microarray.client.dataimport.events.TitleRowChangedEvent;

/**
 * Table to show imported data and color the cells depending on their 
 * status. The cell data can be header, footer, title or actual data.
 * Table also listens for mouse presses for marking column types.
 * 
 * @author mkoski
 *
 */
public class ImportPreviewTable extends JXTable implements MouseMotionListener, MouseListener, 
		ConversionModelChangeListener{
	

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(ImportPreviewTable.class);
	private static final int COLUMN_WIDTH = 120;
	
	/**
	 * Table model for imported data. 
	 * @author mkoski
	 *
	 */
	class ImportPreviewTableModel extends AbstractTableModel {

		private String[] columnTitles;
		private Object[][] arrayData;

		public ImportPreviewTableModel(Object[][] arrayData, String[] columnTitles) {
			this.arrayData = arrayData;
			this.columnTitles = columnTitles;
		}

		public int getRowCount() {
			return arrayData.length;
		}

		public int getColumnCount() {
			return columnTitles.length;
		}

		public String getColumnName(int columnIndex) {
			return columnTitles[columnIndex];
		}

		public Class<?> getColumnClass(int columnIndex) {
			Object firstValue = this.getValueAt(0, columnIndex);
			if (firstValue != null) {
				return firstValue.getClass();
			} else {
				return null;
			}
		}

		/**
		 * Gets value at given row and column. The column count may vary because 
		 * of the headers and footers. If value in given row and column is not set 
		 * the method returns empty string. Otherwise the value of cell is returned
		 * 
		 * @param int row
		 * @param int col
		 * @return value of the cell or empty string if there is no value set
		 */
		public Object getValueAt(int row, int col) {
			if (row >= 0 && row < arrayData.length &&
					col >= 0 && col < arrayData[row].length) {
				return arrayData[row][col];
			} else {
				return "";
			}
		}
	}

	class RowNumberRenderer extends DefaultTableCellRenderer{

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			super.getTableCellRendererComponent(table, value, isSelected, hasFocus,
					row, column);

			// This removes the borders when a cell is focused
			if (hasFocus) {
				Border border = BorderFactory.createEmptyBorder(1,1,1,1);
				setBorder(border);
			}

			this.setBackground(UIManager.getColor("TableHeader.background"));
			this.setHorizontalAlignment(SwingConstants.CENTER);

			return this;
		}
	}


	/**
	 * Renderer which adds posibility to change the color of 
	 * rows which has been marked to header or footer
	 * 
	 * @author mkoski
	 *
	 */
	class ImportCellRenderer extends DefaultTableCellRenderer{

		private int row;
		private int column;
		
		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			super.getTableCellRendererComponent(table, value, isSelected, hasFocus,
					row, column);
			this.row = row;
			this.column = column;
			
			// This removes the borders when a cell is focused
			if (hasFocus) {
				Border border = BorderFactory.createEmptyBorder(1,1,1,1);
				setBorder(border);
			}

			// Highlights
			if(highlightRow && row == mouseOverRow){
				setBackground(new Color(250,250,250));
			} else if(highlightColumn && column == mouseOverColumn){
				setBackground(new Color(250,250,250));
			} else {
				setBackground(UIManager.getColor("Table.background"));
			}

			return this;
		}

		@Override
		protected void paintComponent(Graphics g) {
			
			// Header and footer coloring
			if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
				//Default
				this.setForeground(UIManager.getColor("Table.foreground"));
				
				if (header >= row) {
					this.setForeground(Color.LIGHT_GRAY);
				} 
				
				if(footer <= row){
					this.setForeground(Color.LIGHT_GRAY);
				} 
				
				if(titles == row){
					// This is the column title row
					this.setBackground(Color.LIGHT_GRAY);
					this.setForeground(UIManager.getColor("Table.foreground"));
				} 				
			}
			
			// Do data trimming
			this.setText(screen.getDataTrimmer().doTrimming(this.getText(), column));
			
			// Column coloring on the second step
			if(screen.getCurrentStep() == ImportScreen.Step.SECOND &&	
					columnTypeManager.getColumnCount() > column &&
					columnTypeManager.getColumnType(column) != null 
					&& columnTypeManager.getColumnType(column).getColor() != null){
				
				// Default
				this.setForeground(UIManager.getColor("Table.foreground"));
					
				this.setBackground(columnTypeManager.getColumnType(column).getColor());
				
			}
			
			super.paintComponent(g);
		}
	}
	private static final String HEADER_TOOLTIP = "Click to the row where the header ends";
	private static final String FOOTER_TOOLTIP = "Click to the row where the footer starts";
	private static final String TITLE_TOOLTIP = "Click to the column title row";
	private static final String NO_TOOL_DIALOG_TITLE = "Select a marking tool before clicking the table";
	private static final String NO_TOOL_DIALOG_TEXT = "No tool selected";
	private ConversionModel conversionModel;
	private ImportScreen screen;
	private TableInternalFrame tableFrame;
	private int header;
	private int footer;
	private boolean highlightRow;
	private boolean highlightColumn;
	private int mouseOverRow;
	private int mouseOverColumn;
	private int titles;
	private ColumnTypeManager columnTypeManager;
	private int defaultInitialDelay;
	private int defaultDismissDelay;
	private Point lastDragPoint;

	public ImportPreviewTable(ImportScreen importScreen, TableInternalFrame tableFrame){
		super();
		this.screen = importScreen;
		this.conversionModel = importScreen.getConversionModel();
		this.tableFrame = tableFrame;
		
		this.columnTypeManager = screen.getColumnTypeManager();
		this.setAutoResizeMode(JXTable.AUTO_RESIZE_OFF);
		this.setSortable(false);
		this.setTableHeader(new EditableHeader(this.getColumnModel()));
		this.getTableHeader().setReorderingAllowed(false);
		this.getTableHeader().setEnabled(false);
		
		this.defaultInitialDelay = ToolTipManager.sharedInstance().getInitialDelay();
		this.defaultDismissDelay = ToolTipManager.sharedInstance().getDismissDelay();
		
		// Renderers for data, headers and footers
		setDefaultRenderer(String.class, new ImportCellRenderer());
		setDefaultRenderer(Double.class, new ImportCellRenderer());		

		// Renderer for row numbers
		setDefaultRenderer(Integer.class, new RowNumberRenderer());		

		// Mouse listeners
		this.addMouseListener(this);
		this.addMouseMotionListener(this);
	}
	
	public void createHeaderRenderers(){
		logger.debug("Check the header class of.."+this.getTableHeader()+"..");
		if(!(this.getTableHeader() instanceof EditableHeader)){
			logger.debug("..not editable header, will put a new one");
			this.setTableHeader(new EditableHeader(this.getColumnModel()));
		}
		
		for(int i = 1; i < this.getColumnCount(); i++){
			
			if(!(this.getColumnModel().getColumn(i).getHeaderRenderer() instanceof
					PanelTableHeaderRenderer)){
				
				PanelTableHeaderRenderer renderer = new PanelTableHeaderRenderer(screen,i);

				this.getColumnModel().getColumn(i).setHeaderRenderer(
						renderer);
				
				if(!(this.getColumnModel().getColumn(i) instanceof EditableHeaderTableColumn)){
					((EditableHeader)this.getTableHeader()).recreateTableColumn(columnModel);					
				}

				((EditableHeaderTableColumn)this.getColumnModel().getColumn(i)).
				setHeaderEditor(new PanelTableHeaderEditor(renderer));
				
				this.getColumnModel().getColumn(i).setMinWidth(COLUMN_WIDTH);
			}
			
		}		
		
		screen.getTableFrame().updateAllChipCountComboBoxes();
		screen.getTableFrame().updateAllColumnTitleLabels();	
		this.updateTableHeader();
	}
	
	public PanelTableHeaderRenderer getHeaderRenderer(int columnIndex){
		return (PanelTableHeaderRenderer)this.getColumnModel().getColumn(columnIndex).getHeaderRenderer();
	}

	/**
	 * Sets data to table. This method updates swing components and it made thread safe.
	 * 
	 * @param arrayData
	 * @param columnTitles
	 */
	public void setData(Object[][] arrayData, String[] columnTitles) {
		
		/**
		 * Runnable helper class to set data to table. This is done 
		 * in it's own runnable class because the data must be set 
		 * before column headers can be created. That's why the 
		 * SwingUtilities.invokeAndWait is used because it 
		 * blocks the application until the <code>setModel</code> 
		 * is done.
		 * 
		 * @author mkoski
		 *
		 */
		class SetModelRunnable implements Runnable{

			private Object[][] arrayData;
			private String[] columnTitles;
			
			public SetModelRunnable(Object[][] arrayData, String[] columnTitles) {
				this.arrayData = arrayData;
				this.columnTitles = columnTitles;
			}
			
			public void run() {
				ImportPreviewTable.this.setModel(new ImportPreviewTableModel(arrayData, columnTitles));
				if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
					screen.getTableFrame().addCornerComponents();
				}
			}
			
		}
		
		/**
		 * Makes the rest of table initializing
		 * 
		 * @author mkoski
		 *
		 */
		class InitializeTableRunnable implements Runnable{
			
			public void run(){
				// Take custom renderers in use because previous line sets default 
				// renderers (null), but only for second step 
				if(screen.getCurrentStep() == ImportScreen.Step.SECOND){
					screen.getColumnTypeManager().setColumnCount(getColumnCount());
					createHeaderRenderers();
				}
				
				// Initialize the table's and spinner's footer and header values if needed
				if(screen.getConversionModel().hasColumnTitles() 
						|| screen.getConversionModel().hasFooter()
						|| screen.getConversionModel().hasHeader()){
					// Footer, header or title already set. 
					// Do not initialize the values
				} else {
					initializeHeaderAndFooter();
					tableFrame.resetSpinners();
				}

				// Sets row number column width
				updateRowNumberColumnWidth();
			}
		}
		
		// Set data to table
		try {
			SwingUtilities.invokeAndWait(new SetModelRunnable(arrayData, columnTitles));
		} 
		
		catch (InterruptedException e) { } 
		catch (InvocationTargetException e) { }
		
		// Do other initializations
		SwingUtilities.invokeLater(new InitializeTableRunnable());
	}

	private void updateRowNumberColumnWidth() {	
		if(this.getRowCount()>0 && this.getColumnCount()>0){
			String str = this.getValueAt(this.getRowCount()-1, 0).toString();
			FontMetrics fm = this.getFontMetrics(this.getFont());

			// The number does not fit without a bit extra space (+8)
			int width = SwingUtilities.computeStringWidth(fm, str) + 8;

			this.getColumnModel().getColumn(0).setMaxWidth(width);
			this.getColumnModel().getColumn(0).setMinWidth(width);
			logger.debug("Set row number column width to " + width);
		}
	}

	/**
	 * Updates the cell highlighting according to mouse point
	 * 
	 * @param e Mouse event
	 */
	private void updateHighlight(MouseEvent e){

		int oldMouseOverRow = mouseOverRow;
		int oldMouseOverColumn = mouseOverColumn;

		mouseOverRow = this.rowAtPoint(e.getPoint());
		mouseOverColumn = this.columnAtPoint(e.getPoint());

		// If mouse is off the component, set both to -1 when no highlight
		// is painted
		if(mouseOverRow < 0 || mouseOverColumn < 0){
			mouseOverRow = -1;
			mouseOverColumn = -1;
		}

		// Repainting the table is important to update the highlighted row or column
		// The repainting is a slow operation, so it shouldn't be done everytime 
		// the mouse moves
		if((highlightRow && oldMouseOverRow != mouseOverRow)
				|| (highlightColumn && oldMouseOverColumn != mouseOverColumn)){
			repaint();
		}
	}
	
	/**
	 * Updates table header. Gets the header value from table if cell header row 
	 * is selected. Otherwise the header value is just column numbers.
	 *
	 */
	private void updateTableHeader(){
		TableColumnModelExt model = new DefaultTableColumnModelExt();
		Object[] columnTitles = conversionModel.getColumnTitles();
		
		for(int column = 0; column < columnTitles.length; column++){
			
			if(column >= this.getColumnCount()){
				/*
				 * The TitleRowChangeEvent is fired after data chopping is done. 
				 * This means that this method may be called before the table is updated. 
				 * In this case, if the column count is changed, ArrayIndexOutOfBounds exception 
				 * may occure.
				 *  
				 * So, let's break the loop and wait that this method is called by updateTable 
				 * method a bit later.
				 */
				break;
			}
			TableColumn newColumn = this.getColumnModel().getColumn(column);
			
			if(column >= 1){ 
				
				//Now text to the upper left corner cell
				if(newColumn.getHeaderRenderer() instanceof PanelTableHeaderRenderer){
					//For custom header of second step
					((PanelTableHeaderRenderer)newColumn.getHeaderRenderer()).
					setTitleText(columnTitles[column].toString());
				} else {
					//For step 1				
					newColumn.setHeaderValue(columnTitles[column]);
				}
			}
			
			model.addColumn(newColumn);
		} 
		this.setColumnModel(model);
		
		logger.debug("Table header updated");
	}

	public void setColumnHighlight(boolean enabled){
		highlightColumn = enabled;
		repaint();
	}

	public void setRowHighlight(boolean enabled){
		highlightRow = enabled;
		repaint();
	}

	public void mouseMoved(MouseEvent e) {
		updateHighlight(e);
	}

	public void mouseClicked(MouseEvent e) {
		// Do nothing		
	}

	public void mouseEntered(MouseEvent e) {
		this.setTooltipDelays();
	}

	public void mouseExited(MouseEvent e) {
		updateHighlight(e);
		this.resetTooltipDelays();
	}

	public void mouseDragged(MouseEvent e) {
		if(e.getSource() == this && screen.getCurrentStep() == ImportScreen.Step.FIRST){
			if(tableFrame.isInHeaderMarkingMode()){
				markHeaderAtPoint(e.getPoint());
			} else if(tableFrame.isInFooterMarkingMode()){
				markFooterAtPoint(e.getPoint());
			} else if(tableFrame.isInTitleMarkingMode()){
				markTitleAtPoint(e.getPoint(), false);
			}
		} else {
			markColumnsBetweenPoints(lastDragPoint, e.getPoint());
			if(this.columnAtPoint(e.getPoint()) == -1){
				// User dragged mouse out of the table. Do not set new last point
				return;
			} else {
				lastDragPoint = e.getPoint();
			}
		}
		updateHighlight(e);
	}

	/**
	 * Marks column between two points. This is useful when marking columns 
	 * using mouse drag
	 * 
	 * @param start start point
	 * @param end end point
	 */
	private void markColumnsBetweenPoints(Point start, Point end) {
		int startColumn = this.columnAtPoint(start);
		int endColumn = this.columnAtPoint(end);
		
		if(endColumn < 0){
			// endColumn can get value -1 if user drags mouse out of the table
			endColumn = this.getColumnCount() -2;
		}
		
		int bigger;
		int smaller;
		
		if(startColumn > endColumn){
			bigger = startColumn;
			smaller = endColumn;
		} else {
			bigger = endColumn;
			smaller = startColumn;
		}
		
		for(int column = smaller; column <= bigger; column++){
			markColumn(column);
		}
		
	}

	public void mousePressed(MouseEvent e) {
		// Mark drag start point just in case user starts dragging
		lastDragPoint = e.getPoint();
		
		// Listeners for first step
		if(e.getSource() == this && screen.getCurrentStep() == ImportScreen.Step.FIRST){
			
			// Mark header
			if(tableFrame.isInHeaderMarkingMode()){
				markHeaderAtPoint(e.getPoint());
			} 
			
			// Mark footer
			else if(tableFrame.isInFooterMarkingMode()){
				markFooterAtPoint(e.getPoint());
			} 
			
			// Mark title row
			else if(tableFrame.isInTitleMarkingMode()){
				markTitleAtPoint(e.getPoint(), true);
			} 
			
			else {
				JOptionPane.showMessageDialog(screen.getFrame(), NO_TOOL_DIALOG_TITLE, NO_TOOL_DIALOG_TEXT, JOptionPane.WARNING_MESSAGE);
			}
		} 
		
		// Listeners for second step
		else if(e.getSource() == this && screen.getCurrentStep() == ImportScreen.Step.SECOND){			
			markColumnAtPoint(e.getPoint());
		}
	}
	
	private void markHeaderAtPoint(Point p){
		
		// Set header end row to both table and spinner
		int headerEnd = this.rowAtPoint(p);
		conversionModel.setHeaderEndsRow(headerEnd);
	}
	
	private void markFooterAtPoint(Point p){
		
		// Set footer start row to both table and spinner
		int footerStart = this.rowAtPoint(p);
		conversionModel.setFooterBeginsRow(footerStart);
	}
	
	private void markTitleAtPoint(Point p, boolean removeIfSelected){
		// Sets column title row
		int columnTitleRow = this.rowAtPoint(p);
		
		// Removes title row if the current title row is clicked
		if(titles == columnTitleRow && removeIfSelected){
			conversionModel.setColumnTitleLine(-1);
		} else {
			conversionModel.setColumnTitleLine(columnTitleRow);
		}
		this.repaint();
	}
	
	private void markColumnAtPoint(Point p){
		markColumn(this.columnAtPoint(p));
	}
	
	private void markColumn(int column){
		if(screen.getTableFrame().getSelectedColumnType() != null){
			ColumnType type = screen.getTableFrame().getSelectedColumnType();
			if(column >= 1 && column < screen.getColumnTypeManager().getColumnCount()){					
				
				// Set column type to the column type manager
				screen.getColumnTypeManager().setColumnType(column, type, conversionModel.getCleanColumnTitle(column));
				screen.getColumnTypeManager().setColumnChipNumber(column, screen.getColumnTypeManager().getNextChipNumber(type));
			} else {
				logger.debug("ColumnCount smaller than clicked column index: "+ screen.getColumnTypeManager().getColumnCount());
			}
		}
	}

	public void mouseReleased(MouseEvent e) {
		this.repaint();	
	}

	/**
	 * Sets row number to end the header.
	 * 
	 * @param header Row number to end the header
	 */
	private void setHeaderEndRow(int header){
		this.header = header;
		logger.debug("Header end row set to " + header);
		this.repaint();
	}

	/**
	 * Sets row number to start the footer. Notice that row numbers in the table 
	 * starts from 0, so the value is smaller by one than it is in spinners
	 * 
	 * @param header Row number to start the footer
	 */
	private void setFooterStartRow(int footer){
		if(footer < 0){
			footer = this.getRowCount() +1;
		}
		this.footer = footer;
		logger.debug("Footer start row set to " + footer);
		this.repaint();
	}
	
	private void setColumnTitleRow(int columnTitleRow) {
		this.titles = columnTitleRow;
		this.updateTableHeader();
		this.repaint();
	}

	/**
	 * Resets the header, footer and title rows
	 *
	 */
	public void initializeHeaderAndFooter(){
		// Reset value of this table
		this.header = -1;
		this.footer = this.getRowCount();
		this.titles = -1;
	}
		
	public void decimalSeparatorChanged(DecimalSeparatorChangedEvent e) {
		// Do nothing
	}

	public void delimiterChanged(DelimiterChangedEvent e) {
		// Do nothing
	}

	public void footerChanged(FooterChangedEvent e) {
		this.setFooterStartRow(e.getNewValue());
	}

	public void headerChanged(HeaderChangedEvent e) {
		this.setHeaderEndRow(e.getNewValue());
	}

	public void titleRowChanged(TitleRowChangedEvent e) {
		this.setColumnTitleRow(e.getNewValue());
	}
	
	/**
	 * Makes tooltip to follow the mouse pointer
	 */
	@Override
	public Point getToolTipLocation(MouseEvent event) {
		Point tooltipPoint = new Point(event.getX() + 12, event.getY() + 14);
		return tooltipPoint;
	}
	
	/**
	 * Gets the right tooltip depending on the selected tool and step
	 */
	@Override
	public String getToolTipText() {
		if(screen.getCurrentStep() == ImportScreen.Step.FIRST){
			if(screen.getTableFrame().isInFooterMarkingMode()){
				return FOOTER_TOOLTIP;
			} 
			else if(screen.getTableFrame().isInHeaderMarkingMode()){
				return HEADER_TOOLTIP;
			} 
			else if(screen.getTableFrame().isInTitleMarkingMode()){
				return TITLE_TOOLTIP;
			}
			else {
				return null;
			}
		} else {
			if(screen.getTableFrame().getSelectedColumnType() == null){
				return null;
			} 
			
			else if(screen.getTableFrame().getSelectedColumnType().equals(ColumnType.UNUSED_LABEL)){
				return "Mark as unused";
			}
			
			else {
				return "Select " + screen.getTableFrame().getSelectedColumnType().toString().toLowerCase();
			}
		}
	}
	
	/**
	 * Sets tooltip delays back to default
	 *
	 */
	public void resetTooltipDelays(){
		ToolTipManager.sharedInstance().setInitialDelay(this.defaultInitialDelay);
		ToolTipManager.sharedInstance().setDismissDelay(this.defaultDismissDelay);
	}
	
	/**
	 * Sets tooltips to appear immediately
	 */
	public void setTooltipDelays(){
		ToolTipManager.sharedInstance().setInitialDelay(0);
		ToolTipManager.sharedInstance().setDismissDelay(Integer.MAX_VALUE);
	}

	public void columnTitlesChanged(ColumnTitlesChangedEvent e) {
		// Do nothing
	}

	public void inputFileChanged(InputFileChangedEvent e) {
		// Do nothing
	}
}



