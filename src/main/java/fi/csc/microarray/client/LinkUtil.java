package fi.csc.microarray.client;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.Action;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

public class LinkUtil {
	
	private static final String LINK_WORD = "***";
	
	public static JXHyperlink createLink(String text, Action action) {
		JXHyperlink link = new JXHyperlink();
		link.setBorder(null);
		link.setMargin(new Insets(0, 0, 0, 0));
		link.setAction(action);
		link.setText(text); // must be after setAction
		return link;
	}
	
	public static void addLink(String description, JXHyperlink link, ImageIcon icon, GridBagConstraints c, JComponent component, int maxRowChars, Color background) {
		List<JXHyperlink> list = new LinkedList<JXHyperlink>();
		list.add(link);
		addLinks(description, list, icon, c, component, maxRowChars, background);
	}

	
	public static void addLinks(String description, List<JXHyperlink> links, ImageIcon icon, GridBagConstraints c, JComponent component, int maxRowChars, Color background) {

		String[] words = description.split(" ");
		
		//Split links to separate words to be able to handle links that aren't separated by space character
		List<String> linkSeparatedWords = new LinkedList<String>();
		for (String word : words) {
			if (!word.contains(LINK_WORD)) {
				linkSeparatedWords.add(word + " ");
			} else {
				String wordTail = word + " ";
				
				while (wordTail.length() > 0) {
					int linkWordIndex = wordTail.indexOf(LINK_WORD);
					
					if (linkWordIndex > 0) {
						String wordHead = word.substring(0, linkWordIndex);	
						linkSeparatedWords.add(wordHead);
						wordTail = wordTail.substring(linkWordIndex);
					}
					
					if (linkWordIndex >= 0) {
						String linkWord = wordTail.substring(0, LINK_WORD.length());
						linkSeparatedWords.add(linkWord);
						wordTail = wordTail.substring(LINK_WORD.length());
					}
								
					if (linkWordIndex == -1) {
						linkSeparatedWords.add(wordTail);
						wordTail = "";
					}
				}
			}
		}
				
		int rowChars = 0;

		Iterator<JXHyperlink> linkIterator = links.iterator();
		int rowCount = 0;

		c.gridx = 1;
		c.insets.top = 10;
		JPanel row = null;

		for (int i = 0; i < linkSeparatedWords.size(); i++) {
			if (row == null || rowChars + linkSeparatedWords.get(i).length() > maxRowChars || linkSeparatedWords.get(i).equals("\n ")) {

				FlowLayout flow = new FlowLayout(FlowLayout.LEADING);
				flow.setVgap(0);
				flow.setHgap(0);
				row = new JPanel(flow);
				row.setBackground(background);

				c.gridy++;
				component.add(row, c);
				c.insets.top = 0; // After first row

				rowChars = 0;
				rowCount++;
			}

			if (linkSeparatedWords.get(i).equals(LINK_WORD)) {

				JXHyperlink link = linkIterator.next();
				rowChars += link.getText().length();
				row.add(link);
				
			} else if (!linkSeparatedWords.get(i).equals("\n")) {
				JLabel text = new JLabel(linkSeparatedWords.get(i));
				row.add(text);
				rowChars += linkSeparatedWords.get(i).length();
			}
		}

		c.gridy -= (rowCount - 1);
		c.gridheight = rowCount;
		c.gridx = 0;
		c.insets.top = 15;
		component.add(new JLabel(icon), c);
		c.gridy += (rowCount - 1);
		c.gridheight = 1;
	}

}
