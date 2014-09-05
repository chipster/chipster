package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.LineBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class SessionDetails extends Visualisation implements FocusListener, DocumentListener, MouseListener{
		
	private final String PLEASE_ADD_NOTES = "(Click here to add your notes)";

	private JTextArea notesField;
	private JPanel panel;

	private static final Color BG = Color.white;
	
	public final static int LEFT_WIDTH = 450;
	public final static int INDENTION = 20;

	private JTextArea titleField;

	private JScrollPane scroller;
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	private JPanel getPanelBase(String... layoutArgs) {
		MigLayout layout = null;
		switch (layoutArgs.length) {
		case 0:
			layout = new MigLayout();
			break;			
		case 1:
			layout = new MigLayout(layoutArgs[0]);
			break;			
		case 2:
			layout = new MigLayout(layoutArgs[0], layoutArgs[1]);
			break;
		default:
			layout = new MigLayout(layoutArgs[0], layoutArgs[1], layoutArgs[2]);
			break;
		}
		
		JPanel panel = new JPanel(layout);
		panel.setBackground(BG);
		
		return panel;
	}

	private Component createNotes() {


		// notes visible?
		if (Session.getSession().getPrimaryModule().notesVisibleAtStartup()) {
			
			notesField = new JTextArea();
			notesField.setEditable(true);
			notesField.setLineWrap(true);
			notesField.setWrapStyleWord(true);
			notesField.addFocusListener(this);
			notesField.getDocument().addDocumentListener(this);
			notesField.setRows(1);

			setNotes(application.getSessionNotes());
			setNotesActive(false);
			notesField.setEnabled(true);
			
			notesField.addKeyListener(new KeyListener() {
				
				@Override
				public void keyTyped(KeyEvent e) {
					if (e.getKeyChar() == KeyEvent.VK_ESCAPE) {
						endEditing();
					}
				}
				

				@Override
				public void keyReleased(KeyEvent e) {
				}
				
				@Override
				public void keyPressed(KeyEvent e) {	
				}
			});
			return notesField;
		}
		else {
			return getPanelBase();
		}
	}

	private JTextArea createTitleTextArea(boolean isSingle) {
		
		String name = application.getSessionName();
		
		if (name == null) {
			name = "Unsaved session";
		}
		
		titleField = new JTextArea(name);
		titleField.setFont(titleField.getFont().deriveFont(titleField.getFont().getSize2D() * 1.5f));
		titleField.setEditable(false);
		
		return titleField;			
	}
	
	private void endEditing() {
		scroller.requestFocusInWindow();
	}
	
	private Component createDatasetPanel() {
				
		JPanel panel = getPanelBase("wrap 1, fillx");

		panel.add(createTitleTextArea(true), "growx");			
		panel.add(createNotes(), "gapx " + INDENTION + ", growx");			
								
		return panel;
	}

	@Override
	public JComponent getVisualisation(List<DataBean> datas) throws Exception {
		
		/* See VisualisationTaskManager.VisualisationRunnable.run()
		 * 
		 * This visualization is created so often that we start to see threading
		 * problems. Methods getVisualisation() and visualisationShown() are
		 * synchronized to make things little bit more predictable. 
		 */
		synchronized (this) {

			panel = getPanelBase("gapy 20");						
			scroller = new JScrollPane(panel);			

			return scroller;
		}
	}
	
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean == null;
	}

	@Override
	public boolean canVisualise(List<DataBean> datas) throws MicroarrayException {
		return datas.isEmpty();
	}

	private void setNotes(String text) {

		if (text == null || "".equals(text.trim())) {
			notesField.setText(PLEASE_ADD_NOTES);

		} else {
			notesField.setText(text);			
		}
	}

	private String getNotesContent() {
		if (PLEASE_ADD_NOTES.equals(notesField.getText())) {
			return "";
		} else {
			return notesField.getText();
		}		
	}

	public void focusGained(FocusEvent e) {
		if (e.getComponent() == notesField) {
			setNotesActive(true);
			// user starts writing notes, so remove "please add notes" if needed
			if (PLEASE_ADD_NOTES.equals(notesField.getText())) {
				notesField.setSelectionStart(0);
				notesField.setSelectionEnd(notesField.getText().length());
				notesField.setBorder(new LineBorder(Color.gray));
			}
		}
	}

	public void focusLost(FocusEvent e) {
		// do nothing, content is stored already
		
		if (e.getComponent() == notesField) {
			setNotes(application.getSessionNotes());
			setNotesActive(false);
		}
	}
	
	private void setNotesActive(boolean active) {
		if (active) {
			notesField.setBorder(new LineBorder(Color.gray));
		} else {
			notesField.setBorder(new LineBorder(BG));
			notesField.setEnabled(false);
			notesField.setEnabled(true);
			notesField.setSelectionStart(0);
			notesField.setSelectionEnd(0);
		}
	}	

	public void changedUpdate(DocumentEvent e) {
		fieldUpdated(e);		
	}

	public void insertUpdate(DocumentEvent e) {
		fieldUpdated(e);		
	}

	public void removeUpdate(DocumentEvent e) {
		fieldUpdated(e);		
	}

	private void fieldUpdated(DocumentEvent e) {
		if (e.getDocument() == notesField.getDocument()) {
			setNotesActive(true);
			application.setSessionNotes(getNotesContent());
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		return getVisualisation(new ArrayList<DataBean>());
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		
		endEditing();
	}

	@Override
	public void mouseEntered(MouseEvent e) {
	}

	@Override
	public void mouseExited(MouseEvent e) {
	}

	@Override
	public void mousePressed(MouseEvent e) {
	}

	@Override
	public void mouseReleased(MouseEvent e) {	
	}
	
	@Override
	public void visualisationShown() {
		
		synchronized (this) {

			panel.removeAll();

			// this is EDT, now it's safe to access databeans
			panel.add(createDatasetPanel(), "gapx 20, aligny top, width " + (LEFT_WIDTH + INDENTION));				

			panel.addMouseListener(this);
			scroller.setBorder(null);
		}
	}
}
