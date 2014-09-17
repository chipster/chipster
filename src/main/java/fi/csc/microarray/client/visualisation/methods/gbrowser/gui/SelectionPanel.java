package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import javax.swing.JPanel;
import javax.swing.JTextArea;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;

/**
 * A Simple panel for showing some textual information about the selected items.
 * 
 * @author klemela
 */
public class SelectionPanel extends JPanel implements BrowserSelectionListener {

	private JTextArea textArea;
	private SelectionManager selectionManager;

	public SelectionPanel(SelectionManager selectionManager) {
		
		this.selectionManager = selectionManager;
		selectionManager.addSelectionListener(this);
		this.setLayout(new MigLayout("insets 0", "[grow, fill]", "[grow, fill]"));

		textArea = new JTextArea();
		this.add(textArea);		
	}

	@Override
	public void selectionChanged(DataUrl data, Selectable selectable, Object source) {
		textArea.setText(selectionManager.getSelectionText());
	}
}
