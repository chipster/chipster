package fi.csc.microarray.client.dataimport.table;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

import fi.csc.microarray.client.dataimport.ImportScreen;


/**
 * Custom renderer for the header of the JTable, which contains two JLabels
 * and a comboBox to show and edit column information. Placing interactive
 * components in the header needs lot of tweaking which can be done with 
 * classes EditableHeader, EditableHeaderTableColumn and EditableHeaderUI.
 * Just showing those components with this panel isn't enough, they need also 
 * an editor class, which is implemented in PanelTableHeaderEditor.
 * 
 * @author klemela
 */
public class PanelTableHeaderRenderer  extends JPanel 
	implements TableCellRenderer, ItemListener  {

	private JComboBox chipCombo;
	private JLabel titleLabel;
	private JLabel typeLabel;
	
	private ImportScreen screen;
	private int columnIndex;

	public PanelTableHeaderRenderer(ImportScreen screen, int i) {
		this.screen = screen;
		columnIndex = i;
		
		chipCombo = new JComboBox();
		titleLabel = new JLabel("Title");
		typeLabel = new JLabel("Type");
		
		chipCombo.addItemListener(this);
		//titleLabel.setForeground(UIManager.getColor("Label.disabledForeground"));
		
		
		this.setLayout(new BorderLayout());	
		
		Border border = UIManager.getBorder("TableHeader.cellBorder");
		border.getBorderInsets(this).set(5,5,5,5);
		setBorder(border);
				
		setBackground(UIManager.getColor("TableHeader.background"));

		add(titleLabel,BorderLayout.NORTH);
		add(typeLabel,BorderLayout.WEST);
		add(chipCombo,BorderLayout.EAST);
	}

	public void setSelectedIndex(int index) {
		//logger.debug("Modelsize: " + chipCombo.getModel().getSize()+" index: " +index);
		if(chipCombo.getModel().getSize()>index && index >= 0){
			chipCombo.setSelectedIndex(index);
		}
	}

	public int getSelectedIndex() {
		return chipCombo.getSelectedIndex();    
	}

	public JComboBox getCombo(){
		return chipCombo;
	}
	
	public void setTitleText(String text){
		titleLabel.setText(text);
		this.update();
	}
	
	public void setTypeText(String text){
		typeLabel.setText("<html><b>" + text + "</b></html>");
		this.update();
	}
	
	public void setTypeColor(Color c){
		typeLabel.setForeground(c);
		this.update();
	}
	
	public void update(){		
		//repaint on this panel isn't enough, this will do litle bit more
		screen.getTableFrame().getTable().getTableHeader().repaint();
	}
	
	public void setTypeToolTipText(String text){
		typeLabel.setToolTipText(text);
	}

	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {
		if (value instanceof Integer) {
			setSelectedIndex(((Integer)value).intValue());
		}    
		return this;
	}
	
	public void itemStateChanged(ItemEvent e) {
		if(e.getSource() instanceof JComboBox){
			Integer newValue = (Integer)e.getItem();
			
			// Change the chip number only if it has really changed.
			if(newValue != screen.getColumnTypeManager().getColumnChipNumber(columnIndex)){
				screen.getColumnTypeManager().setColumnChipNumber(columnIndex, newValue);
			}
		}
	}
}

