package fi.csc.microarray.client.dataview;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.DataItem;

/**
 * A GUI component for showing and editing text notes about each dataset.
 * It has two JTextArea components, one non-editable (intended for "final"
 * kind of data, such as the name of the operation which produced the
 * current dataset, along with all operation parameters) and one editable
 * (for the user's personal notes). Both are placed together in one
 * JScrollPane, so that the length of the notes is not limited.
 * 
 * @author Janne KÃ¤ki
 *
 */
public class DetailsPanel extends JPanel implements PropertyChangeListener, FocusListener, ActionListener, DocumentListener {

	public static final int MINIMISED_HEIGHT = 30;
	private final Color NOTES_BACKGROUND = Color.WHITE;
	private final String PLEASE_ADD_NOTES = "Add your notes here...";
	private final String NO_DATASET_SELECTED = "No dataset selected";

	private final String CL_MINIMISED = "minimised-panel";
	private final int  CL_MINIMISED_INDEX = 0;
	private final String CL_MAIN = "main-panel";
	private final int CL_MAIN_INDEX = 1;
	
	private JTextArea attributesField = new JTextArea();  // not editable by user
	private JLabel nameField = new JLabel();  // not editable by user
	private JTextArea notesField = new JTextArea();
	private JLabel minimisedNotesField = new JLabel();// not editable by user
	
	private JButton hideButton = new JButton("Hide");
	private JButton showButton = new JButton("Show");
	private JScrollPane detailsScroller;
	private JPanel panelChanger = new JPanel(new CardLayout());
	private Dataset currentData;
	private ClientApplication application = Session.getSession().getApplication();
	
	JComponent parent;
	
	/**
	 * Initializes a new DetailsPanel.
	 */
	public DetailsPanel(JComponent parent) {
		this.parent = parent;

		// create static part of details field
		nameField.setForeground(VisualConstants.DETAILS_NAME_FOREGROUND_COLOR);
		attributesField.setEditable(false);
		attributesField.setLineWrap(true);
		attributesField.setWrapStyleWord(true);
		attributesField.setForeground(VisualConstants.DETAILS_ATTRIBUTES_FOREGROUND_COLOR);
		attributesField.setBackground(VisualConstants.TEXTAREA_UNEDITABLE_BACKGROUND);
		minimisedNotesField.addFocusListener(this);
		
		// created dynamic (user editable) part
		notesField.setEditable(true);
		notesField.setLineWrap(true);
		notesField.setWrapStyleWord(true);
		notesField.setBackground(NOTES_BACKGROUND);		
		notesField.addFocusListener(this);
		notesField.getDocument().addDocumentListener(this);
		
		// make it scrollable
        JScrollPane notesScroller = new JScrollPane(notesField);
        notesScroller.setBorder(BorderFactory.createEmptyBorder());
        notesScroller.setPreferredSize(new Dimension(200, 20));

        // make main panel
        JPanel contentPane = new JPanel(new BorderLayout());
        JPanel upperPanel = new JPanel(new BorderLayout());
        hideButton.addActionListener(this);
        JPanel upperUpperPanel = new JPanel(new BorderLayout());
        upperUpperPanel.add(hideButton, BorderLayout.EAST);
        upperUpperPanel.add(nameField, BorderLayout.CENTER);
        upperPanel.add(upperUpperPanel, BorderLayout.NORTH);
        upperPanel.add(attributesField, BorderLayout.CENTER);
        contentPane.add(upperPanel, BorderLayout.NORTH);
        contentPane.add(notesScroller, BorderLayout.CENTER);
        contentPane.setPreferredSize(new Dimension(240, VisualConstants.DETAILS_PANEL_HEIGHT - 20));

        // make scroller panel (main)
		detailsScroller = new JScrollPane(contentPane);
		detailsScroller.setBorder(null);
        detailsScroller.setHorizontalScrollBarPolicy(
                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		detailsScroller.setVerticalScrollBarPolicy(
				ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
		detailsScroller.setPreferredSize(new Dimension(
                VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.DETAILS_PANEL_HEIGHT));

		// make minimised panel
		JPanel minimisedPanel = new JPanel(new BorderLayout());			
		JPanel innerMinimisedPanel = new JPanel(new BorderLayout());
		innerMinimisedPanel.setMinimumSize(new Dimension(0,0));
		showButton.addActionListener(this);
		innerMinimisedPanel.add(showButton, BorderLayout.EAST);
		innerMinimisedPanel.add(minimisedNotesField, BorderLayout.WEST);
		minimisedPanel.add(innerMinimisedPanel, BorderLayout.NORTH);
		minimisedPanel.setPreferredSize(new Dimension(
				VisualConstants.LEFT_PANEL_WIDTH, MINIMISED_HEIGHT));
		minimisedPanel.setMinimumSize(new Dimension(0,0));
		
		// make panel changer (minimised/main)
		panelChanger.add(minimisedPanel, CL_MINIMISED); // ORDER OF ADD(...)'S IS SIGNIFICANT
		panelChanger.add(detailsScroller, CL_MAIN); // ORDER OF ADD(...)'S IS SIGNIFICANT
		
		// put it all together
		this.setLayout(new BorderLayout());
		this.add(panelChanger, BorderLayout.CENTER);
		this.setPreferredSize(panelChanger.getComponent(0).getPreferredSize());
		this.setMinimumSize(new Dimension(0,0));
		this.disable();
		
		// start listening
		application.addPropertyChangeListener(this);
	}
	
	private String getNameText() {
		if (currentData != null) {
			return currentData.getOperation().getDefinition().getFullName();
			
		} else {			
			return null;
		}		
	}
	
	/**
	 * @return A String containing descriptions of the chosen dataset's
	 * 		   attributes - that is, name, date, and details about the
	 * 		   operation (including parameters) that produced it.
	 */
	private String getAttributeText() {
		if (currentData != null) {
			StringBuffer attrib = new StringBuffer();
			attrib.append(currentData.getDate().toString());
			Operation lastOper = currentData.getOperation();
			if (lastOper != null) {
				
				attrib.append("\nOperation: " + lastOper.getCategoryName() +
						" / " +lastOper.getID());
				
				List<Parameter> params = lastOper.getParameters();
				if (params != null) {
					attrib.append("\n");
					for (int i = 0; i < params.size(); i++) {
						attrib.append(params.get(i).toString());
						if (i != params.size()-1) {
							attrib.append(", ");
						}
					}
				}
			}
			return attrib.toString();
			
		} else {
			return null;
		}
	}
    
    public void setViewedData(DataItem data) {
    	
        if (data != null && data instanceof Dataset) {
            Dataset dataBean = (Dataset)data;
            this.currentData = dataBean;
            nameField.setText(getNameText());
            attributesField.setText(getAttributeText());
            setNotes(dataBean.getNotes());
            notesField.setEnabled(true);
            notesField.setBackground(NOTES_BACKGROUND);
            
        } else {
            disable();
        }
    }
	
	public void disable() {
		this.currentData = null;
		minimisedNotesField.setText(NO_DATASET_SELECTED);
		nameField.setText(NO_DATASET_SELECTED);
		attributesField.setText("");
		notesField.setText("");		
		notesField.setEnabled(false);
		notesField.setBorder(null);
	}
	
	private void setNotes(String text) {
		minimisedNotesField.setText(this.getNameText());
		
		if (text == null || "".equals(text.trim())) {
			notesField.setText(PLEASE_ADD_NOTES);
			
		} else {
			notesField.setText(text);			
		}
	}
	
	private String getNotes() {
		if (PLEASE_ADD_NOTES.equals(notesField.getText())) {
			return "";
		} else {
			return notesField.getText();
		}		
	}
    
    private void showPanel(String name, int index) {
        ((CardLayout)panelChanger.getLayout()).show(panelChanger, name);
        this.setPreferredSize(panelChanger.getComponent(index).getPreferredSize());
        parent.validate();
    }
    
	public void focusGained(FocusEvent e) {
		// user starts writing notes, so remove "please add notes" if needed
		if (PLEASE_ADD_NOTES.equals(notesField.getText())) {
			notesField.setText("");
		}		
	}

	public void focusLost(FocusEvent e) {
		// do nothing, content is already stored
	}
    
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == hideButton) {
			showPanel(CL_MINIMISED, CL_MINIMISED_INDEX);					
		} else if (e.getSource() == showButton) {
			showPanel(CL_MAIN, CL_MAIN_INDEX);
		}
	}
    
    public void propertyChange(PropertyChangeEvent dataEvent) {
        DataItem data = application.getSelectionManager().getSelectedItem();        
        this.setViewedData(data);
    }

	public void changedUpdate(DocumentEvent e) {
		notesFieldUpdated();		
	}

	public void insertUpdate(DocumentEvent e) {
		notesFieldUpdated();		
	}

	public void removeUpdate(DocumentEvent e) {
		notesFieldUpdated();		
	}
	
	private void notesFieldUpdated() {
		if (this.currentData != null) {
			this.currentData.setNotes(getNotes()); // update notes
		}
	}
}
