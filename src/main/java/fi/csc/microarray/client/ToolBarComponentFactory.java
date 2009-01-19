package fi.csc.microarray.client;

import java.awt.Color;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JSpinner;
import javax.swing.JToggleButton;

public class ToolBarComponentFactory {
	
	public static void initialiseButton(AbstractButton button, boolean leftBorder, boolean rightBorder){
		button.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createMatteBorder(0,leftBorder?1:0,0,rightBorder?1:0,Color.LIGHT_GRAY),
				BorderFactory.createEmptyBorder(4,4,4,4)));
	}
	
	//J B U T T O N S //////////////////////////////////
	
	public static JButton createButton(String text, ImageIcon icon, boolean leftBorder, boolean rightBorder){
		JButton button = new JButton(text,icon);		
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JButton createButton(String text, boolean leftBorder, boolean rightBorder){
		JButton button = new JButton(text);
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JButton createButton(ImageIcon icon, boolean leftBorder, boolean rightBorder){
		JButton button = new JButton(icon);
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JButton createButton(boolean leftBorder, boolean rightBorder){
		JButton button = new JButton();
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	// J T O G G L E  B U T T O N S //////////////////////////////////////////////////
	
	public static JToggleButton createToggleButton(String text, ImageIcon icon, boolean leftBorder, boolean rightBorder){
		JToggleButton button = new JToggleButton(text,icon);
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JToggleButton createToggleButton(String text, boolean leftBorder, boolean rightBorder){
		JToggleButton button = new JToggleButton(text);
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JToggleButton createToggleButton(ImageIcon icon, boolean leftBorder, boolean rightBorder){
		JToggleButton button = new JToggleButton(icon);
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	public static JToggleButton createToggleButton(boolean leftBorder, boolean rightBorder){
		JToggleButton button = new JToggleButton();
		initialiseButton(button,leftBorder,rightBorder);
		return button;
	}
	
	// J S P I N N E R //////////////////////////////////////////////
	
	public static JSpinner createSpinner(){
		JSpinner spinner = new JSpinner();
		spinner.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createMatteBorder(0,1,0,1,Color.LIGHT_GRAY),
				BorderFactory.createEmptyBorder(0,0,0,0)));
		return spinner;
	}
	
	// J C O M B O  B O X /////////////////////////////////////////////
	
	public static JComboBox createComboBox() {
		JComboBox combo = new JComboBox();
		combo.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createMatteBorder(0,0,0,0,Color.LIGHT_GRAY),
				BorderFactory.createEmptyBorder(0,0,0,0)));
		return combo;
	}
	
	// J C H E C K  B O X /////////////////////////////////////////////
	
	public static JCheckBox createCheckBox(String text){
		JCheckBox cbox = new JCheckBox(text);
		cbox.setOpaque(false);
		return cbox;
	}
}
