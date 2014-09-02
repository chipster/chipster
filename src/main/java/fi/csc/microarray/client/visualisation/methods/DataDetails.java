package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Collection;
import java.util.List;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.LineBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import net.miginfocom.swing.MigLayout;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.ClientFocusTraversalPolicy;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationToolBar;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class DataDetails extends Visualisation implements FocusListener, DocumentListener, MouseListener{
	
	class LinkMouseListener implements MouseListener {
		private JXHyperlink link;

		public LinkMouseListener(JXHyperlink link) {
			this.link = link;
		}

		@Override
		public void mouseClicked(MouseEvent e) {
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			link.setOpaque(true);
		}

		@Override
		public void mouseExited(MouseEvent e) {
			link.setOpaque(false);
		}

		@Override
		public void mousePressed(MouseEvent e) {
		}

		@Override
		public void mouseReleased(MouseEvent e) {
		}				
	}

	private final String PLEASE_ADD_NOTES = "(Click here to add your notes)";

	private JTextArea notesField;
	private JPanel panel;

	private static final Color BG = Color.white;

	public static final String COMMAND = "command";
	public static final String RENAME_COMMAND = "rename";
	
	public final static int LEFT_WIDTH = 450;
	public final static int INDENTION = 20;

	private List<DataBean> datas;

	private JTextArea titleField;

	private JScrollPane scroller;

	private Vector<JComponent> focusableLinks = new Vector<JComponent>();

	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	private JPanel getPanelBase() {
		return getPanelBase(null);
	}
	private JPanel getPanelBase(String layoutArgs) {
		JPanel panel = new JPanel(new MigLayout(layoutArgs));
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
			//notesField.setColumns(LEFT_WIDTH / notesField.getFont().getSize()); //not accurate
			notesField.setRows(1);

			setNotes(datas.get(0).getNotes());
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

	private Component createVisualisations() {
		JPanel panel = getPanelBase("wrap 1");
		
		List<VisualisationMethod> visualisations = VisualisationToolBar.getMethodsFor(datas);
			
		for (VisualisationMethod method : visualisations) {
	
			JXHyperlink link = new JXHyperlink();
			
			link.setBackground(new Color(0.95f, 0.95f, 0.95f));					
			
			link.addMouseListener(new LinkMouseListener(link));
			//hide focus border because it doesn't obey component size
			link.setFocusPainted(false);

			link.setIcon(method.getIcon());
			link.setIconTextGap(16);
			link.addActionListener(new VisualisationStarter(method, Session.getSession().getApplication()));
			link.setText(method.getName());
			
			panel.add(link, "width 300!, height 50!");
			
			focusableLinks.add(link);
		}
		
		return panel;
	}

	private JTextArea createTitleTextArea(DataBean data, boolean isSingle) {
		
		titleField = new JTextArea(data.getName());
		titleField.setFont(titleField.getFont().deriveFont(titleField.getFont().getSize2D() * 1.5f));

		// rename functionality
		if (isSingle) {
			titleField.addFocusListener(this);

			titleField.getDocument().addDocumentListener(this);
			setTitleActive(false);
			
			titleField.addKeyListener(new KeyListener() {
				
				@Override
				public void keyTyped(KeyEvent e) {
					if (e.getKeyChar() == KeyEvent.VK_ENTER || e.getKeyChar() == KeyEvent.VK_ESCAPE) {
						titleField.setText(titleField.getText().replace("\n", ""));
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

		} else {
			titleField.setEditable(false);
		}
		
		return titleField;			
	}
	
	private void endEditing() {
		scroller.requestFocusInWindow();
	}
	
	private Component createDatasetPanel(List<DataBean> datas) {
				
		JPanel panel = getPanelBase("wrap 1");
		
		if (datas.size() == 1) {

			DataBean data = datas.get(0);

			panel.add(createTitleTextArea(data, true), "growx");			
			panel.add(createDateLabel(data), "gapx " + INDENTION);
			panel.add(createNotes(), "gapx " + INDENTION + ", growx");			
			panel.add(createToolLabel(data), ", gapy 20");
			createParameterTable(panel);
			
		} else {

			for (DataBean data : datas) {
				panel.add(createTitleTextArea(data, false));				
				panel.add(createToolLabel(data));							
			}
		}
						
		return panel;
	}

	private JLabel createDateLabel(DataBean data) {
		JLabel dateLabel = new JLabel(data.getDate().toString());
		return dateLabel;		
	}

	private Component createToolLabel(DataBean data) {
		JLabel toolLabel = new JLabel(getNameText(data));
		toolLabel.setFont(toolLabel.getFont().deriveFont(Font.BOLD));
		return toolLabel;
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

			// this is not in EDT, do only minimal stuff to avoid problems if something else modifies databeans at the same time
			this.datas = datas;

			panel = getPanelBase("gapy 20");						
			scroller = new JScrollPane(panel);			

			return scroller;
		}
	}
	
	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean != null;
	}

	@Override
	public boolean canVisualise(List<DataBean> datas) throws MicroarrayException {
		return datas.size() > 0;
	}

	private String getNameText(DataBean data) {
		if (data != null) {
			return data.getOperationRecord().getFullName();

		} else {			
			return null;
		}		
	}

	/**
	 * @param panel2 
	 * @return A String containing descriptions of the chosen dataset's
	 * 		   attributes - that is, name, date, and details about the
	 * 		   operation (including parameters) that produced it.
	 */
	private JPanel createParameterTable(JPanel panel) {

		final int TOOL_WIDTH = 250;
		
			
		OperationRecord operationRecord = datas.get(0).getOperationRecord();			
		if (operationRecord != null) {

			Collection<ParameterRecord> params = operationRecord.getParameters();

			if (params != null) {
				for (ParameterRecord parameterRecord : params) {

					// find out default value and human readable value
					OperationDefinition tool = application.getOperationDefinition(operationRecord.getNameID().getID());
					
					String defaultValue = null;
					String valueString = null;
										
					if (tool != null) {
						defaultValue = tool.getParameterDefaultValue(parameterRecord);
						valueString = tool.getHumanReadableParameterValue(parameterRecord);
						
					} else {						
						valueString = parameterRecord.getValue();
					}

					JTextArea name = new JTextArea(parameterRecord.getNameID().getDisplayName());					
					JTextArea value =  new JTextArea(valueString);
					
					name.setLineWrap(true);
					value.setLineWrap(true);
					
					name.setWrapStyleWord(true);
					value.setWrapStyleWord(true);
					
					name.setEditable(false);
					value.setEditable(false);

					// fade out default values
					if (defaultValue != null && defaultValue.equals(parameterRecord.getValue())) {
						value.setForeground(Color.gray);
					}
					
					panel.add(name, "split 2, gapx " + INDENTION + ", width " + TOOL_WIDTH);
					panel.add(value, "growx, pushx, wrap");
				}
			}
		}

		return panel;
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
		} else if (e.getComponent() == titleField) {
			setTitleActive(true);
		}
	}

	public void focusLost(FocusEvent e) {
		// do nothing, content is stored already
		
		if (e.getComponent() == notesField) {
			setNotes(datas.get(0).getNotes());
			setNotesActive(false);
 
		} else if (e.getComponent() == titleField) {
			setTitleActive(false);
			
			if ("".equals(titleField.getText().trim())) {
				titleField.setText(datas.get(0).getName());
			}
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
	
	private void setTitleActive(boolean active) {
		if (active) {
			titleField.setBorder(new LineBorder(Color.gray));
		} else {
			titleField.setBorder(new LineBorder(BG));
			titleField.setEnabled(false);
			titleField.setEnabled(true);
			titleField.setSelectionStart(0);
			titleField.setSelectionEnd(0);
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
			if (this.datas != null) {
				this.datas.get(0).setNotes(getNotesContent()); // update notes
			}
		} else if (e.getDocument() == titleField.getDocument()) {
			setTitleActive(true);
			if (!"".equals(titleField.getText().trim())) {
				application.renameDataItem(datas.get(0), titleField.getText().trim());
			}
		}
	}
	
	@Override
	public boolean isForMultipleDatas() {
		return true;
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		// not used because this is for multiple datas
		return null;
	}
	
	private class VisualisationStarter implements ActionListener {
		
		private VisualisationMethod method;
		private ClientApplication application;

		private VisualisationStarter(VisualisationMethod method, ClientApplication clientApplication) {
			this.method = method;
			this.application = clientApplication;
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			application.setVisualisationMethod(method, null, datas, getFrame().getType());
		}
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
	
	public void updateFocusTraversal(JComponent base){
		
		Vector<Component> order = new Vector<Component>();
		order.addAll(focusableLinks);
		
		base.addFocusListener(new FocusListener() {
			
			@Override
			public void focusLost(FocusEvent arg0) {
			}
			
			@Override
			public void focusGained(FocusEvent arg0) {
				if (focusableLinks.size() > 0) {
					focusableLinks.get(0).requestFocusInWindow();
				}
			}
		});
	
		base.setFocusCycleRoot(true);
		base.setFocusTraversalPolicy(new ClientFocusTraversalPolicy(order));
	}
	
	@Override
	public void visualisationShown() {
		
		synchronized (this) {

			panel.removeAll();

			// this is EDT, now it's safe to access databeans
			panel.add(createDatasetPanel(datas), "gapx 20, aligny top, width " + (LEFT_WIDTH + INDENTION));				
			panel.add(createVisualisations(), "aligny top");

			panel.addMouseListener(this);
			scroller.setBorder(null);
			updateFocusTraversal(scroller);

			if (getFrame().getVariables() != null) {
				for (Variable variable : getFrame().getVariables()) {
					//These are set by rename menu command
					if (COMMAND.equals(variable.getName()) && RENAME_COMMAND.equals(variable.getExpression())) {
						titleField.requestFocusInWindow();
						titleField.getCaret().setDot(titleField.getText().length());
					}
				}
			}
		}
	}
}
