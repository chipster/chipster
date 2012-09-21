package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
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
import fi.csc.microarray.client.visualisation.VisualisationMethodRepository.VisualisationMethodOrderComparator;
import fi.csc.microarray.client.visualisation.VisualisationToolBar;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule.VisualisationMethods;

public class DataDetails extends Visualisation implements FocusListener, DocumentListener, MouseListener{

	private final String PLEASE_ADD_NOTES = "(Add your notes here)";

	private JTextArea notesField;
	private JPanel panel;

	private static final Color BG = Color.white;

	public static final String COMMAND = "command";
	public static final String RENAME_COMMAND = "rename";
	
	final int LEFT_WIDTH = 400;
	final int INDENTION = 20;

	private List<DataBean> datas;

	private JPanel cachePanel;

	private JTextArea titleField;

	private JScrollPane scroller;

	private Vector<JComponent> focusableLinks = new Vector<JComponent>();

	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	private JPanel getPanelBase(String title, boolean indent) {
		JPanel panel = new JPanel(new BorderLayout());
		panel.setBackground(BG); 
		if (title != null) {
			
			JLabel titleLabel = new JLabel(" " + title);
			//titleLabel.setForeground(Color.gray);
			titleLabel.setFont(titleLabel.getFont().deriveFont(Font.BOLD));
			panel.add(titleLabel, BorderLayout.NORTH);
			//panel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder(caption));
		}
		if (indent) {
			JPanel indention = new JPanel();
			indention.setBackground(BG);
			indention.setPreferredSize(new Dimension(INDENTION, 1));
			panel.add(indention, BorderLayout.WEST);
		}
		return panel;
	}
	

	private JPanel getPanelBase(String title) {
		return getPanelBase(title, true);
	}

	private Component createNotes() {


		// notes visible?
		if (Session.getSession().getPrimaryModule().notesVisibleAtStartup()) {

			JPanel notesPanel = getPanelBase(null);
			
			JPanel panel = new JPanel(new BorderLayout());
			panel.setBackground(BG);
			
			notesField = new JTextArea();
			notesField.setEditable(true);
			notesField.setLineWrap(true);
			notesField.setWrapStyleWord(true);
			notesField.addFocusListener(this);
			notesField.getDocument().addDocumentListener(this);
			notesField.setColumns(LEFT_WIDTH / notesField.getFont().getSize()); //not accurate
			notesField.setRows(1);
			panel.add(notesField, BorderLayout.CENTER);
			panel.add(new JLabel(" "), BorderLayout.SOUTH);
			notesPanel.add(panel, BorderLayout.CENTER);

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
			return notesPanel;
		}
		else {
			return getPanelBase(null);
		}

	}

	private Component createVisualisations() {
		JPanel visualisationsPanel = getPanelBase(" ");

		JPanel panel = new JPanel(new GridBagLayout());
		panel.setBackground(BG);
		GridBagConstraints c = new GridBagConstraints();
		
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.WEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.ipadx = 10;
		c.insets = new Insets(5, 5, 5, 5);
		
		List<VisualisationMethod> visualisations = VisualisationToolBar.getMethodsFor(datas);
		visualisations.remove(VisualisationMethod.NONE);
		visualisations.remove(VisualisationMethods.DATA_DETAILS);
		
		LinkedList<VisualisationMethod> orderedMethods = new LinkedList<VisualisationMethod>(visualisations);
		Collections.sort(orderedMethods, new VisualisationMethodOrderComparator());
		
		for (VisualisationMethod method : orderedMethods) {
			
//			JLabel icon = new JLabel(resizeImage(method.getIcon()));
			JLabel icon = new JLabel(method.getIcon());
			JXHyperlink link = new JXHyperlink();
			link.addActionListener(new VisualisationStarter(method, Session.getSession().getApplication()));
			link.setText(method.getName());
			c.gridx = 0;
			c.weightx = 0;
			panel.add(icon, c);
			c.gridx = 1;
			c.weightx = 1.0;
			panel.add(link, c);
			c.gridy++;
			
			focusableLinks.add(link);
		}
		
		visualisationsPanel.add(panel, BorderLayout.CENTER);
		
		return visualisationsPanel;
	}

	private Component createParameters(DataBean data, boolean isSingle) {

		JPanel parametersPanel = getPanelBase(getNameText(data));
		
		if (isSingle) {
			parametersPanel.add(getParameterTable(), BorderLayout.CENTER);
		}

		return parametersPanel;
	}

	private Component createDatasetDetails(DataBean data, boolean isSingle) {
		
		JPanel datasetPanel = getPanelBase(null);

		titleField = new JTextArea(data.getName());
		titleField.setFont(titleField.getFont().deriveFont(titleField.getFont().getSize2D() * 1.5f));

		if (isSingle) {
			titleField.addFocusListener(this);
			//titleField.getDocument().putProperty("filterNewlines", Boolean.TRUE);
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
		
		datasetPanel.add(titleField, BorderLayout.NORTH);
		
		if (isSingle) {
			JPanel panel = new JPanel(new GridBagLayout());
			panel.setBackground(BG);

			JLabel dateLabel = new JLabel(data.getDate().toString());
			JLabel locationLabel = new JLabel("Location: chipster.csc.fi ");

			cachePanel = new JPanel(new BorderLayout());
			cachePanel.setBackground(BG);


			JLabel cacheStartLabel = new JLabel("(");
			cachePanel.add(cacheStartLabel, BorderLayout.WEST);

			JXHyperlink link = new JXHyperlink();
			link.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					cachePanel.removeAll();
					cachePanel.add(new JLabel("(cached locally)"));
					cachePanel.validate();
				}
			});
			link.setText("Get local copy");
			cachePanel.add(link, BorderLayout.CENTER);

			JLabel cacheEndLabel = new JLabel(")");
			cachePanel.add(cacheEndLabel, BorderLayout.EAST);

			//"Local file"
			//"Remote file (chipster.csc.fi)"
			//"Remote file (chipster.csc.fi) with local copy"

			GridBagConstraints c = new GridBagConstraints();
			c.gridx = 0;
			c.gridy = 0;
			c.gridwidth = 3;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.weightx = 0;
			c.weighty = 0;

			panel.add(dateLabel, c);

			c.gridy++;
			c.gridx = 0;
			c.gridwidth = 1;
			panel.add(locationLabel, c);

			c.gridx++;
			panel.add(cachePanel, c);

			c.gridx++;
			c.weightx = 1.0;
			JPanel xSpaceFiller = new JPanel();
			xSpaceFiller.setBackground(BG);
			panel.add(xSpaceFiller, c);
			c.gridx = 0;
			c.gridy++;
			c.weighty = 1.0;
			c.weightx = 0;
			JPanel ySpaceFiller = new JPanel();
			ySpaceFiller.setBackground(BG);
			panel.add(ySpaceFiller, c);

			datasetPanel.setPreferredSize(new Dimension(LEFT_WIDTH + INDENTION, locationLabel.getFont().getSize() * 5));
			datasetPanel.add(panel, BorderLayout.CENTER);
		}
		return datasetPanel;
	}
	
	private void endEditing() {
		scroller.requestFocusInWindow();
	}
	
	private Component createDatasetPanel(List<DataBean> datas) {
		
		JPanel panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		panel.setBackground(BG);
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.NONE;
		c.ipady = 10;

		if (datas.size() == 1) {

			DataBean data = datas.get(0);

			panel.add(createDatasetDetails(data, true), c);
			c.gridy++;

			panel.add(emptyIfMultipleDatas(createNotes()), c);
			c.gridy++;

			panel.add(emptyIfMultipleDatas(createParameters(data, true)), c);
			c.gridy++;
			
		} else {

			for (DataBean data : datas) {
				Component datasetDetails = createDatasetDetails(data, false);
				panel.add(datasetDetails, c);
				c.gridy++;
				
				panel.add(createParameters(data, false), c);
				c.gridy++;
				
				panel.setPreferredSize(new Dimension(
						INDENTION + LEFT_WIDTH, 
						//multiply by 3 leave some space for the tool name
						(int) (datasetDetails.getPreferredSize().getHeight() * 3 + panel.getPreferredSize().getHeight())));
			}
		}
		
		c.gridy++;
		addEmptyRow(panel, c, -1);
		
		c.gridy = 0;
		c.gridx++;
		addEmptyColumn(panel, c, -1);
						
		return panel;
	}

	@Override
	public JComponent getVisualisation(List<DataBean> datas) throws Exception {
		this.datas = datas;

		panel = new JPanel();
		panel.setLayout(new GridBagLayout());
		panel.setBackground(BG);
		
		final int TOP_MARGIN = 10;
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.NONE;
		c.ipady = 10;
		
		addEmptyColumn(panel, c, 20);

		c.gridx++;
		
		addEmptyRow(panel, c, TOP_MARGIN);
		
		c.gridy++;
		panel.add(createDatasetPanel(datas), c);
		
		c.gridy++;
		addEmptyRow(panel, c, -1);
		
		c.gridy = 0;
		c.gridx++;
		addEmptyColumn(panel, c, 50);

		c.gridy = 0;
		c.gridx++;
		//addEmptyRow(panel, c, TOP_MARGIN);
		
		c.gridy++;
		c.gridheight = 3;
		panel.add(createVisualisations(), c);
		c.gridy += c.gridheight;
		c.gridheight = 1;
		
		addEmptyRow(panel, c, -1);
		
		c.gridy = 0;
		c.gridx++;
		addEmptyColumn(panel, c, -1);
		
		panel.addMouseListener(this);
		
		scroller = new JScrollPane(panel);
		updateFocusTraversal(scroller);

		return scroller;
	}
	
	private Component emptyIfMultipleDatas(Component c) {
		if (datas.size() == 1) {
			return c;
		} else {
			JPanel empty = new JPanel();
			empty.setBackground(BG);
			empty.setPreferredSize(new Dimension(LEFT_WIDTH + INDENTION, 1));
			return empty;
		}
	}
		
	private void addEmptyColumn(JPanel container, GridBagConstraints c, int width) {
		JPanel panel = new JPanel();
		panel.setBackground(BG);
		if (width > 0) {
			panel.setPreferredSize(new Dimension(width, 1));
			container.add(panel, c);
		} else {
			double tmp_weightx = c.weightx;
			c.weightx = 1.0;
			container.add(panel, c);
			c.weightx = tmp_weightx;
		}
	}
	
	private void addEmptyRow(JPanel container, GridBagConstraints c, int height) {
		JPanel panel = new JPanel();
		panel.setBackground(BG);
		if (height > 0) {
			panel.setPreferredSize(new Dimension(1, height));
			container.add(panel, c);
		} else {
			double tmp_weighty = c.weighty;
			c.weighty = 1.0;
			container.add(panel, c);
			c.weighty = tmp_weighty;
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
	 * @return A String containing descriptions of the chosen dataset's
	 * 		   attributes - that is, name, date, and details about the
	 * 		   operation (including parameters) that produced it.
	 */
	private JPanel getParameterTable() {
		JPanel panel = new JPanel();
		panel.setBackground(BG);
		panel.setLayout(new GridBagLayout());
		
		final int TOOL_WIDTH = 300;
		final int VALUE_WIDTH = LEFT_WIDTH - TOOL_WIDTH;
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1.0;
		c.weighty = 0;
		
		if (datas != null) {
			
			OperationRecord operationRecord = datas.get(0).getOperationRecord();			
			if (operationRecord != null) {
				
				Collection<ParameterRecord> params = operationRecord.getParameters();
				
				if (params != null) {
					for (ParameterRecord parameterRecord : params) {

						OperationDefinition tool = application.getOperationDefinition(operationRecord.getNameID().getID());
						String defaultValue = tool.getParameter(parameterRecord.getNameID().getID()).getValueAsString();
						
						JLabel name = new JLabel(parameterRecord.getNameID().getDisplayName());
						JLabel value = new JLabel(parameterRecord.getValue());
						
						if (defaultValue.equals(parameterRecord.getValue())) {
							value.setForeground(Color.gray);
						}
						
						int height = name.getFont().getSize() + 8;
						
						name.setPreferredSize(new Dimension(TOOL_WIDTH, height));
						value.setPreferredSize(new Dimension(VALUE_WIDTH, height));
						
						c.gridx = 0;
						panel.add(name, c);
						c.gridx = 1;
						panel.add(value, c);
						
						c.gridy++;
					}
				}
			}
			c.weighty = 1.0;
			
			JPanel spaceFiller = new JPanel();
			spaceFiller.setBackground(BG);
			
			panel.add(spaceFiller, c);
			
			return panel;
			//return attrib.toString();

		} else {
			return null;
		}
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
				//notesField.setText("");
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
		// TODO Auto-generated method stub
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
