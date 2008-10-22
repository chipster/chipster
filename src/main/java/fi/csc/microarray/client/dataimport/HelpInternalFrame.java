package fi.csc.microarray.client.dataimport;

import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JLabel;

import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

/**
 * Help panel for Import tool.
 * 
 * @author pklemela
 *
 */
public class HelpInternalFrame extends SimpleInternalFrame {
	
	private JLabel textContent = new JLabel();
	
	public static final int HELP_HEIGHT = 190;	
	
	/**
	 * Help text for the first step
	 */
	private static final String TEXT_FIRST_STEP =
"<html><body><table><tr>"+
"<th>1</th><th>2</th><th>3</th><th>4</th><th>5</th>"+
"</tr><tr>"+
"<td>Select the delimiter to split the data table in columns. Empty or unused columns can be ignored later.</td>"+
"<td>Change decimal separator if needed. Decimal separator can't be the same as the column delimiter character.</td>"+
"<td>Select the tool Mark header, and click the last row of the header (the row just before the data).</td>"+
"<td>Select the tool Mark footer, and click the first row of footer (the row just after the data).</td>"+
"<td>If the data table has a row for column titles, it can be copied to the table header by clicking the row with tool Mark title row.</td>"+
"</tr></table></body></html>";
	
/**
 * Help text for the second step
 */
private static final String TEXT_SECOND_STEP =
"<html><body><table><tr>"+
"<th>1</th><th>2</th><th>3</th><th>4</th><th>5</th>"+
"</tr><tr>"+
"<td>Select the tool sample from the toolbar and mark sample columns by clicking them with the mouse.</td>"+
"<td>If the file contains other data colums, mark them with corresponding toolbar buttons. A pattern can be copied to the whole table with guess -button.</td>"+
"<td>Check that all data columns from the same chip have the same chip number in the column title.</td>"+
"<td>Chip counts -panel shows the count of selected columns of each type.</td>"+
"<td>If the flag column isn't already in P/M/A format, change the values to corresponding characters with flag modification.</td>"+
"</tr></table></body></html>";

	
	public HelpInternalFrame(){
		super("Help");
		textContent.setMinimumSize(new Dimension(0,0));
		this.setContent(textContent);
		this.setBackground(Color.WHITE);
	}
	
	/**
	 * Shows help text for current step.
	 * 
	 * @param step
	 */
	public void showHelp(ImportScreen.Step step){
		if(step == ImportScreen.Step.FIRST){
			textContent.setText(TEXT_FIRST_STEP);
		} else {
			textContent.setText(TEXT_SECOND_STEP);
		}		
	}
}
