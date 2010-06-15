package fi.csc.microarray.client;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.dialog.RenameDialog;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.Dataset.Link;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;

/**
 * Context menu for set of selected DataItems.
 * 
 * @see fi.csc.microarray.databeans.DataItem
 * 
 * @author Mikko Koski, Petri KlemelÃ¤, Aleksi Kallio
 */
public class ClientContextMenu extends JPopupMenu implements ActionListener, PopupMenuListener {

	private static final Logger logger = Logger.getLogger(ClientContextMenu.class);

	private static class LinkInfo {
		private Dataset source;
		private Dataset target;
		private Link link;

		public LinkInfo(Dataset source, Dataset target, Link link) {
			this.source = source;
			this.target = target;
			this.link = link;
		}

		public Dataset getSource() {
			return source;
		}

		public Dataset getTarget() {
			return target;
		}

		public Link getLink() {
			return link;
		}
	}

	private SwingClientApplication application;

	private DataItem selectedItem;

	private JMenuItem visualiseMenuItem;
	private JMenuItem renameMenuItem;
	private JMenuItem deleteMenuItem;
	private JMenuItem importMenuItem;
	private JMenuItem exportMenuItem;
	private JMenuItem historyMenuItem;
	private JMenuItem saveWorkflowItem;
	private JMenu linksMenu;
	private JMenu linkToMenu;
	private JMenu unlinkMenu;
	private JMenu metadataLinkMenu;

	/**
	 * Possible destinations of link which may be created
	 */

	/**
	 * Maps linkDestinations and data beans
	 */
	private Map<JMenuItem, LinkInfo> linkMap = new HashMap<JMenuItem, LinkInfo>();
	private Map<JMenuItem, LinkInfo> unlinkMap = new HashMap<JMenuItem, LinkInfo>();

	public ClientContextMenu(SwingClientApplication application) {
		this.application = application;

		this.addPopupMenuListener(this);

		visualiseMenuItem = new JMenuItem("Visualise");
		visualiseMenuItem.setFont(this.getFont().deriveFont(Font.BOLD));
		metadataLinkMenu = new JMenu("Link to phenodata");
		metadataLinkMenu.setIcon(VisualConstants.LINK_PHENODATA_MENUICON);
		linksMenu = new JMenu("Links between selected");
		linkToMenu = new JMenu("Link");
		linkToMenu.setIcon(VisualConstants.LINK_MENUICON);
		linksMenu.add(linkToMenu);
		unlinkMenu = new JMenu("Unlink");
		unlinkMenu.setIcon(VisualConstants.UNLINK_MENUICON);
		linksMenu.add(unlinkMenu);

		renameMenuItem = new JMenuItem("Rename");
		deleteMenuItem = new JMenuItem("Delete");
		deleteMenuItem.setIcon(VisualConstants.DELETE_MENUICON);
		importMenuItem = new JMenuItem("Import files...");
		exportMenuItem = new JMenuItem("Export...");
		exportMenuItem.setIcon(VisualConstants.EXPORT_MENUICON);
		historyMenuItem = new JMenuItem("View history as text");
		historyMenuItem.setIcon(VisualConstants.GENERATE_HISTORY_ICON);
		saveWorkflowItem = new JMenuItem("Save workflow");
		
		
		visualiseMenuItem.addActionListener(this);
		renameMenuItem.addActionListener(this);
		deleteMenuItem.addActionListener(this);
		importMenuItem.addActionListener(this);
		exportMenuItem.addActionListener(this);
		historyMenuItem.addActionListener(this);
		saveWorkflowItem.addActionListener(this);

	}

	public void setOptionsFor(List<DataItem> items) {

		// To avoid wrong menu when the selections is only null
		while (items.contains(null)) {
			items.remove(null);
		}

		boolean multipleSelected = items.size() > 1;
		if (items != null && items.size() > 0) {
			selectedItem = items.get(0); // This is data for single selection
		}
		this.removeAll();

		// Items for the datasets
		if (items != null && items.size() > 0) {
			this.add(visualiseMenuItem);
			this.addSeparator();
			this.add(metadataLinkMenu);
			this.add(linksMenu);
			this.addSeparator();
			this.add(renameMenuItem);
			this.add(deleteMenuItem);
			this.addSeparator();
			this.add(saveWorkflowItem);
			this.addSeparator();
			this.add(exportMenuItem);
			this.add(historyMenuItem);

			if (multipleSelected) {
				// multiple selected DataBeans
				this.visualiseMenuItem.setEnabled(application.isSelectedDataVisualisable());

				this.renameMenuItem.setEnabled(false);
				this.exportMenuItem.setEnabled(true);
				this.historyMenuItem.setEnabled(false);
				this.saveWorkflowItem.setEnabled(false);

			} else {
				if (selectedItem instanceof Dataset) {
					// single selected DataBean
					this.visualiseMenuItem.setEnabled(true);
					this.renameMenuItem.setEnabled(true);
					this.exportMenuItem.setEnabled(true);
					this.historyMenuItem.setEnabled(true);
					
					// workflow items only enabled for normalised data
					boolean normalisedDataSelected = ChipsterInputTypes.GENE_EXPRS.isTypeOf((Dataset)selectedItem);
					this.saveWorkflowItem.setEnabled(normalisedDataSelected);
				} else {
					// single selected DataFolder
					this.visualiseMenuItem.setEnabled(false);
					this.renameMenuItem.setEnabled(true);
					this.exportMenuItem.setEnabled(true);
					this.historyMenuItem.setEnabled(false);
					this.saveWorkflowItem.setEnabled(false);
				}
			}

		} else {
			// Items for empty area
			this.add(importMenuItem);
		}

		// update link options
		updateLinkMenus(items);
	}

	private void updateLinkMenus(List<DataItem> items) {
		
		// do we have exactly two beans (normal linking)
		if (items.size() == 2 && items.get(0) != null && items.get(0) instanceof Dataset && items.get(1) != null && items.get(1) instanceof Dataset) {

			this.metadataLinkMenu.setEnabled(false);
			this.linksMenu.setEnabled(true);

			Dataset one = (Dataset) items.get(0);
			Dataset other = (Dataset) items.get(1);

			linkToMenu.removeAll();

			createNewLinkMenu(one, other);
			createRemoveLinkMenu(one, other);

		// do we have exactly one bean (phenodata shortcut linking)
		} else if (items.size() == 1 && items.get(0) != null && items.get(0) instanceof Dataset) {

			this.metadataLinkMenu.setEnabled(true);
			this.linksMenu.setEnabled(false);

			int itemsCreated = 0;
			
			// phenodata cannot link to phenodata
			if (!((Dataset) items.get(0)).queryFeatures("/phenodata").exists()) {

				List<Dataset> allDatas = application.getAllDataBeans();

				allDatas.remove(selectedItem); // cannot link to itself

				// gather annotating beans
				for (Dataset data : allDatas) {
					if (data.queryFeatures("/phenodata").exists()) {				
						JMenuItem annotateMenuItem = new JMenuItem(getBeanText(data));
						metadataLinkMenu.add(annotateMenuItem);
						linkMap.put(annotateMenuItem, new LinkInfo((Dataset)selectedItem, data, Link.ANNOTATION));
						annotateMenuItem.addActionListener(this);
						itemsCreated++;
					}
				}
			} 
			
			// if we found nothing to annotate with, disable menu
			if (itemsCreated == 0) {
				this.metadataLinkMenu.setEnabled(false);
			}
			
		// set of selected items is not linkable
		} else {
			this.metadataLinkMenu.setEnabled(false);
			this.linksMenu.setEnabled(false);
		}
	}

	private void createRemoveLinkMenu(Dataset one, Dataset other) {

		this.unlinkMenu.setEnabled(false);

		for (LinkInfo link : getSharedLinks(one, other)) {
			this.unlinkMenu.setEnabled(true);

			JMenuItem unlink = new JMenuItem(link.link.toString() + ": " + getLinkText(one, other, false));
			unlinkMap.put(unlink, link);
			this.unlinkMenu.add(unlink);
			unlink.addActionListener(this);
		}
	}

	private String getBeanText(Dataset one) {
		return one.getName();
	}

	private void createNewLinkMenu(Dataset one, Dataset other) {
		Dataset.Link[] links = Dataset.Link.userEditableValues();

		for (Dataset.Link link : links) {
			JMenu linkMenu = new JMenu(link.toString());
			linkToMenu.add(linkMenu);

			JMenuItem oneMenuItem = new JMenuItem(getLinkText(one, other, false));
			linkMenu.add(oneMenuItem);
			linkMap.put(oneMenuItem, new LinkInfo(one, other, link));
			oneMenuItem.addActionListener(this);

			JMenuItem otherMenuItem = new JMenuItem(getLinkText(one, other, true));
			linkMenu.add(otherMenuItem);
			linkMap.put(otherMenuItem, new LinkInfo(other, one, link));
			otherMenuItem.addActionListener(this);
		}
	}

	private String getLinkText(Dataset one, Dataset other, boolean reversed) {
		//Unicode characters didn't work in windows
		//String arrow = reversed ? " \u2190 " : " \u2192 ";
		String arrow = !reversed ? "  ->  " : "  <-  ";
		return getBeanText(one) + arrow + getBeanText(other);				
	}

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == visualiseMenuItem) {
			application.visualiseWithBestMethod(FrameType.MAIN);

		} else if (source == renameMenuItem) {
			new RenameDialog(application, this.selectedItem);

		} else if (source == deleteMenuItem) {
			if (selectedItem instanceof DataFolder) {
				application.deleteDatas((DataFolder) selectedItem);
			} else {
				application.deleteDatas(application.getSelectionManager().getSelectedDataBeans().toArray(new DataItem[0]));
			}
		} else if (source == importMenuItem) {
			try {
				application.openFileImport();
			} catch (Exception me) {
				application.reportException(me);
			}

		} else if (source == exportMenuItem) {
			application.exportSelectedItems();
		} else if (source == historyMenuItem) {
			// Don't do for the folders
			if (selectedItem instanceof Dataset) {
				logger.debug("History activated");
				Dataset data = (Dataset) this.selectedItem;
				if (data != null) {
					application.showHistoryScreenFor(data);
				}
			}
		} else if (source == saveWorkflowItem) {
			application.saveWorkflow();
		}
		
		
		else if (linkMap.keySet().contains(source)) {
			// Create link
			// Note! The actual link direction is opposite to graph
			// arrow directions
			LinkInfo linkData = linkMap.get(source);
			application.createLink(linkData.getTarget(), linkData.getSource(), linkData.getLink());

		} else if (unlinkMap.keySet().contains(source)) {
			// Remove link
			LinkInfo connector = unlinkMap.get((JMenuItem) source);

			logger.debug("Removing link from " + connector.getSource());
			printLinksForData((Dataset) selectedItem);

			application.removeLink(connector.getSource(), connector.getTarget(), connector.getLink());
		}
	}

	public void popupMenuCanceled(PopupMenuEvent e) {
		// Nothing to do

	}

	public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {
		// Nothing to do
	}

	public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
		if (selectedItem instanceof Dataset) {
			historyMenuItem.setEnabled(true);
		} else {
			historyMenuItem.setEnabled(false);
		}
	}

	private static List<LinkInfo> getSharedLinks(Dataset one, Dataset other) {

		List<LinkInfo> sharedLinks = new ArrayList<LinkInfo>();

		for (Link linkType : Link.values()) {

			for (Dataset target : one.getLinkTargets(linkType)) {
				if (other == target) {
					sharedLinks.add(new LinkInfo(one, target, linkType));
				}
			}

			for (Dataset target : other.getLinkTargets(linkType)) {
				if (one == target) {
					sharedLinks.add(new LinkInfo(other, target, linkType));
				}
			}
		}

		// Sort
		Collections.sort(sharedLinks, new Comparator<LinkInfo>() {
			public int compare(LinkInfo o1, LinkInfo o2) {
				return o1.getLink().toString().compareTo(o2.getLink().toString());
			}

		});
		return sharedLinks;
	}

	public static void printLinksForData(Dataset data) {
		for (Link type : Link.values()) {
			for (Dataset source : data.getLinkSources(type)) {
				logger.debug(source.getName() + " <" + type + "> " + data.getName());
			}
		}
		for (Link type : Link.values()) {
			for (Dataset target : data.getLinkTargets(type)) {
				logger.debug(data.getName() + " <" + type + "> " + target.getName());
			}
		}
	}
}
