package fi.csc.microarray.client.waiting;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseMotionAdapter;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.border.EmptyBorder;

/**
 * Makes UI look disabled and consume all user input.
 * 
 * @author Aleksi Kallio
 *
 */
public class WaitGlassPane extends JComponent {
	private static final int BG_ALPHA = 150;
	private static final int MSG_ALPHA = 200;
	
	private JLabel messageLabel = new JLabel();

	public WaitGlassPane() {
		
		// make glass pane white and use alpha channel
		setOpaque(false);
		Color background = new Color(Color.WHITE.getRed(), Color.WHITE.getGreen(), Color.WHITE.getBlue(), BG_ALPHA);
		setBackground(background);
		
		// add message label
		setLayout(new GridBagLayout());
		add(messageLabel, new GridBagConstraints());
		messageLabel.setOpaque(true);
		messageLabel.setBorder(new EmptyBorder(10, 10, 10, 10));

		// make this consume all user input
		addMouseListener(new MouseAdapter() {});
		addMouseMotionListener(new MouseMotionAdapter() {});
		addKeyListener(new KeyListener() {
			public void keyPressed(KeyEvent e) {
				e.consume();
			}
			public void keyReleased(KeyEvent e) {
				e.consume();
			}
			public void keyTyped(KeyEvent e) {
				// ignore
			}			
		});
		setFocusTraversalKeysEnabled(false);
	}

	@Override
	protected void paintComponent(Graphics g) {
		g.setColor(getBackground());
		g.fillRect(0, 0, getSize().width, getSize().height);
	}

	@Override
	public void setBackground(Color background) {
		super.setBackground(background);
		Color messageBackground = new Color(background.getRed(), background.getGreen(), background.getBlue(), MSG_ALPHA);
		messageLabel.setBackground(messageBackground);
	}

	public void startWaiting(String message) {
		messageLabel.setVisible(true);
		messageLabel.setText(message);
		messageLabel.setForeground(getForeground());
		setVisible(true);
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		requestFocusInWindow();
	}

	public void stopWaiting() {
		setCursor(null);
		setVisible(false);
	}
}
