package fi.csc.microarray.util;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.LineBorder;

import fi.csc.microarray.constants.VisualConstants;

/**
 * Mimicks JDK 1.6 SplashScreen-class.
 * 
 * @author akallio
 */
public class SplashScreen {
	
	private JFrame frame;
	private String[] lines = new String[] { "", "", "" , ""};
	private int TEXT_HEIGHT = 70;
	private JLabel textLabel;
	private JPanel textPanel;
	
	public SplashScreen(Icon icon) {
		frame = new JFrame();
		frame.setLayout(null);
		frame.setUndecorated(true);		
		frame.setSize(icon.getIconWidth(), icon.getIconHeight()+TEXT_HEIGHT);
		frame.setLocationRelativeTo(null);  
		frame.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		frame.getRootPane().setBorder(new LineBorder(VisualConstants.SPLASH_BORDER_COLOR, 1));
		
		JLabel imageLabel = new JLabel(icon);
		imageLabel.setSize(icon.getIconWidth(), icon.getIconHeight());
		frame.add(imageLabel);
		imageLabel.setLocation(0, 0);
		textLabel = new JLabel();
		textLabel.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
		
		textLabel.setFont(VisualConstants.SPLASH_SCREEN_FONT);
		
		textPanel = new JPanel(new BorderLayout());
		textPanel.setLocation(0, icon.getIconHeight());
		textPanel.setSize(icon.getIconWidth(), TEXT_HEIGHT);
		textPanel.setBackground(Color.WHITE);
		textPanel.setOpaque(true);
		textPanel.setBorder(BorderFactory.createMatteBorder(1, 0, 0, 0, VisualConstants.SPLASH_BORDER_COLOR));
		textPanel.add(textLabel);
		
		frame.add(textPanel);
		frame.setVisible(true);
	}

	public void close() {
		frame.setVisible(false);
		frame.dispose();
	}

	public void write(String string) {
		lines[lines.length-1] = lines[lines.length-1] + string; 
		updateText();
	}
	
	public void writeLine(String newLine) {
		for (int i = 0; i < lines.length-1; i++) {
			lines[i] = lines[i+1];
		}
		
		lines[lines.length-1] = newLine;
		
		updateText();
	}

	private void updateText() {
		String s = "<html>";
		for (String line : lines) {
			s += (line + "<br>");
		}
		textLabel.setText(s + "</html>");
	}
}
