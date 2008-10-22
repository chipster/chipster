package fi.csc.microarray.util;

import javax.swing.JTextPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.HyperlinkEvent.EventType;

public class BrowsableHtmlPanel {
	
	public static JTextPane createHtmlPanel() {
		JTextPane htmlPanel = new JTextPane(); 
		htmlPanel.setEditable(false);
		htmlPanel.setContentType("text/html");
		htmlPanel.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent e) {
				if (e.getEventType() == EventType.ACTIVATED) {
					try {
						BrowserLauncher.openURL(e.getURL().toString());            
					} catch (Exception ioe) {
						ioe.printStackTrace();
					}
				}
			}
		});
		
		return htmlPanel;
	}

}
