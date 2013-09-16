package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;

public interface BrowserSelectionListener {

	void selectionChanged(DataUrl data, Selectable selectable, Object source);
}
