package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;

/**
 * A Simple panel for showing some textual information about the selected items.
 * 
 * @author klemela
 */
public class SelectionPanel extends JPanel implements BrowserSelectionListener {

	private JPanel panel;
	private JTextArea textArea;
	private SelectionManager selectionManager;

	public SelectionPanel(SelectionManager selectionManager) {
		
		this.selectionManager = selectionManager;
		selectionManager.addSelectionListener(this);
		
		panel = new JPanel();
		panel.setLayout(new BorderLayout());

		textArea = new JTextArea();
		panel.add(textArea, BorderLayout.CENTER);

		JScrollPane scroller = new JScrollPane(panel);
		scroller.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		this.setLayout(new BorderLayout());
		this.add(scroller, BorderLayout.CENTER);			
	}

	@Override
	public void selectionChanged(DataUrl data, Selectable selectable, Object source) {
		textArea.setText(selectionManager.getSelectionText());
	}
}
