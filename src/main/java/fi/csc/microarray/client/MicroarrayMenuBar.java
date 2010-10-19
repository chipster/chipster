package fi.csc.microarray.client;

import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.URL;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.dialog.ClipboardImportDialog;
import fi.csc.microarray.client.dialog.RenameDialog;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationToolBar;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.util.Files;

@SuppressWarnings("serial")
public class MicroarrayMenuBar extends JMenuBar implements PropertyChangeListener {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(MicroarrayMenuBar.class);

	private SwingClientApplication application;

	private JMenu fileMenu = null;
	private JMenu importMenu = null;
	private JMenuItem directImportMenuItem = null;
	private JMenuItem importFromURLMenuItem = null;
	private JMenuItem importFromClipboardMenuItem = null;
	private JMenuItem openWorkflowsMenuItem = null;
	private JMenuItem addDirMenuItem = null;
	private JMenuItem exportMenuItem = null;
	private JMenuItem quitMenuItem = null;
	private JMenu editMenu = null;
	private JMenuItem renameMenuItem = null;
	private JMenuItem deleteMenuItem = null;
	private JMenu viewMenu = null;
	private JMenuItem restoreViewMenuItem = null;
	private JMenu fontSizeMenu = null;
	private JMenu workflowsMenu = null;
	private JMenu helpInfoMenu = null;
	private JMenuItem aboutMenuItem = null;
	private JMenuItem contentMenuItem;
	private JMenuItem startedMenuItem;
	private JMenuItem saveWorkflowMenuItem;
	private JMenuItem helpWorkflowMenuItem;
	private JMenuItem saveSnapshotMenuItem;
	private JMenu recentWorkflowMenu;
	private JMenuItem loadSnapshotMenuItem;
	private JMenuItem taskListMenuItem;
	private JMenuItem clearSessionMenuItem;
	private JMenuItem selectAllMenuItem;
	private JMenuItem historyMenuItem;
	private JMenuItem maximiseVisualisationMenuItem;
	private JMenuItem visualiseMenuItem;
	private JMenuItem detachMenuItem;
	private JMenu visualisationMenu;
	private JMenu openRepoWorkflowsMenu;

	private boolean hasRepoWorkflows;

	public MicroarrayMenuBar(SwingClientApplication application) {
		this.application = application;
		add(getFileMenu());
		add(getEditMenu());
		add(getViewMenu());
		add(getWorkflowsMenu());
		add(getHelpInfoMenu());
		application.addPropertyChangeListener(this);
	}

	public void updateMenuStatus() {

		logger.debug("updating menubar when selected is " + application.getSelectionManager().getSelectedItem());
		DataSelectionManager selectionManager = application.getSelectionManager();

		DataBean selectedDataBean = selectionManager.getSelectedDataBean();
		boolean somethingSelected = selectionManager.getSelectedItem() != null;
		boolean workflowCompatibleDataSelected = false;

		if (selectedDataBean != null) {
			workflowCompatibleDataSelected = Session.getSession().getPrimaryModule().isWorkflowCompatible(selectedDataBean);
		}

		historyMenuItem.setEnabled(selectedDataBean != null && application.getSelectionManager().getSelectedDataBeans().size() == 1);

		visualiseMenuItem.setEnabled(selectedDataBean != null);

		VisualisationMethod method = application.getVisualisationFrameManager().getFrame(FrameType.MAIN).getMethod();
		visualisationMenu.setEnabled(method != null && method != VisualisationMethod.NONE);

		recentWorkflowMenu.setEnabled(workflowCompatibleDataSelected);
		openWorkflowsMenuItem.setEnabled(workflowCompatibleDataSelected);
		openRepoWorkflowsMenu.setEnabled(workflowCompatibleDataSelected && hasRepoWorkflows);

		saveWorkflowMenuItem.setEnabled(workflowCompatibleDataSelected);

		exportMenuItem.setEnabled(somethingSelected);
		renameMenuItem.setEnabled(somethingSelected);
		deleteMenuItem.setEnabled(somethingSelected);
	}

	/**
	 * This method initializes fileMenu
	 * 
	 * @return javax.swing.JMenu
	 */
	private JMenu getFileMenu() {
		if (fileMenu == null) {
			fileMenu = new JMenu();
			fileMenu.setText("File");
			fileMenu.setMnemonic('F');
			fileMenu.add(getDirectImportMenuItem());
			fileMenu.add(getAddDirMenuItem());
			fileMenu.add(getImportMenu());
			fileMenu.addSeparator();
			fileMenu.add(getExportMenuItem());
			fileMenu.addSeparator();
			fileMenu.add(getLoadSnapshotMenuItem());
			fileMenu.add(getSaveSnapshotMenuItem());
			fileMenu.add(getClearSessionMenuItem());
			fileMenu.addSeparator();
			fileMenu.add(getQuitMenuItem());
		}
		return fileMenu;
	}

	private JMenu getImportMenu() {
		if (importMenu == null) {
			importMenu = new JMenu();
			importMenu.setText("Import from");
			
			if (!application.isStandalone) {
				Module primaryModule = Session.getSession().getPrimaryModule();
				primaryModule.addImportMenuItems(importMenu);
			}
				
			importMenu.add(getImportFromURLMenuItem());
			importMenu.add(getImportFromClipboardMenuItem());
		}
		return importMenu;
	}

	private JMenuItem getDirectImportMenuItem() {
		if (directImportMenuItem == null) {
			directImportMenuItem = new JMenuItem();
			directImportMenuItem.setText("Import files...");
			directImportMenuItem.setAccelerator(KeyStroke.getKeyStroke('I', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			directImportMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						application.openFileImport();
					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return directImportMenuItem;
	}

	private JMenuItem getImportFromURLMenuItem() {
		if (importFromURLMenuItem == null) {
			importFromURLMenuItem = new JMenuItem();
			importFromURLMenuItem.setText("URL...");
			importFromURLMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						application.openURLImport();
					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return importFromURLMenuItem;
	}

	private JMenuItem getHelpWorkflowMenuItem() {
		if (helpWorkflowMenuItem == null) {
			helpWorkflowMenuItem = new JMenuItem();
			helpWorkflowMenuItem.setText("More information...");
			helpWorkflowMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						application.viewHelp("chipster-manual/workflows.html#ready");
						
					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return helpWorkflowMenuItem;
	}
	
	private JMenuItem getOpenWorkflowMenuItem() {
		if (openWorkflowsMenuItem == null) {
			openWorkflowsMenuItem = new JMenuItem();
			openWorkflowsMenuItem.setText("Run...");
			openWorkflowsMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						File workflow = application.openWorkflow();
						addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));

					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return openWorkflowsMenuItem;
	}

	private JMenu getOpenRepositoryWorkflowMenu() {
		if (openRepoWorkflowsMenu == null) {

			// Load repository content
			String[][] repoWorkflows = Session.getSession().getPrimaryModule().getRepositoryWorkflows();
			this.hasRepoWorkflows = repoWorkflows.length > 0;

			openRepoWorkflowsMenu = new JMenu();
			openRepoWorkflowsMenu.setText("Run from Chipster repository");

			// Add repository content
			for (String[] flow : repoWorkflows) {
				JMenuItem item = new JMenuItem(flow[0]);
				item.addActionListener(new RepoWorkflowActionListener(flow[0], flow[1]));
				openRepoWorkflowsMenu.add(item);
			}

			// Add help
			openRepoWorkflowsMenu.addSeparator();
			openRepoWorkflowsMenu.add(getHelpWorkflowMenuItem());
			
		}
		return openRepoWorkflowsMenu;
	}

	private class RepoWorkflowActionListener implements ActionListener {

		private String resourceName;
		private String name;

		public RepoWorkflowActionListener(String name, String resourceName) {
			this.name = name;
			this.resourceName = resourceName;
		}

		public void actionPerformed(ActionEvent e) {
			URL url = getClass().getResource(resourceName);
			((SwingClientApplication) application).runWorkflow(url);
			addRecentWorkflow(name, url);
		}
	}

	public void addRecentWorkflow(String name, URL url) {
		if (url != null) { // if the fileChooser is cancelled
			// Check if this exists already
			for (int i = 0; i < recentWorkflowMenu.getItemCount(); i++) {
				JMenuItem menuItem = recentWorkflowMenu.getItem(i);
				if (menuItem.getText().equals(name)) {
					recentWorkflowMenu.remove(menuItem);
				}

			}
			recentWorkflowMenu.add(createRunWorkflowMenuItem(name, url));
		}
	}

	private JMenuItem getImportFromClipboardMenuItem() {
		if (importFromClipboardMenuItem == null) {
			importFromClipboardMenuItem = new JMenuItem();
			importFromClipboardMenuItem.setText("Clipboard...");
			importFromClipboardMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						new ClipboardImportDialog(Session.getSession().getApplication());
					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return importFromClipboardMenuItem;
	}

	private JMenuItem getExportMenuItem() {
		if (exportMenuItem == null) {
			exportMenuItem = new JMenuItem();
			exportMenuItem.setText("Export dataset(s) or folder...");
			exportMenuItem.setAccelerator(KeyStroke.getKeyStroke('E', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			exportMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.exportSelectedItems();
				}
			});
		}
		return exportMenuItem;
	}

	private JMenuItem getQuitMenuItem() {
		if (quitMenuItem == null) {
			quitMenuItem = new JMenuItem();
			quitMenuItem.setText("Quit");
			quitMenuItem.setAccelerator(KeyStroke.getKeyStroke('Q', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			quitMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.quit();
				}
			});
		}
		return quitMenuItem;
	}

	private JMenu getEditMenu() {
		if (editMenu == null) {
			editMenu = new JMenu();
			editMenu.setText("Edit");
			editMenu.setMnemonic('E');
			editMenu.add(getRenameMenuItem());
			editMenu.add(getDeleteMenuItem());
			editMenu.add(getHistoryMenuItem());
			editMenu.addSeparator();
			editMenu.add(getSelectAllMenuItem());
		}
		return editMenu;
	}

	private JMenuItem getHistoryMenuItem() {
		if (historyMenuItem == null) {
			historyMenuItem = new JMenuItem();
			historyMenuItem.setText("Show history...");
			historyMenuItem.setIcon(VisualConstants.GENERATE_HISTORY_ICON);
			historyMenuItem.setAccelerator(KeyStroke.getKeyStroke('H', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			historyMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					DataBean data = application.getSelectionManager().getSelectedDataBean();
					if (data != null) {
						application.showHistoryScreenFor(data);
					}
				}
			});
		}
		return historyMenuItem;
	}

	private JMenuItem getRenameMenuItem() {
		if (renameMenuItem == null) {
			renameMenuItem = new JMenuItem();
			renameMenuItem.setText("Rename selected item...");
			renameMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F2, 0));
			renameMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					new RenameDialog(application, application.getSelectionManager().getSelectedItem());
				}
			});
		}
		return renameMenuItem;
	}

	private JMenuItem getDeleteMenuItem() {
		if (deleteMenuItem == null) {
			deleteMenuItem = new JMenuItem();
			deleteMenuItem.setText("Delete selected item");
			deleteMenuItem.setIcon(VisualConstants.DELETE_MENUICON);
			deleteMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0));
			deleteMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.deleteDatas(application.getSelectionManager().getSelectedDataBeans().toArray(new DataItem[0]));
				}
			});
		}
		return deleteMenuItem;
	}

	private JMenuItem getSelectAllMenuItem() {
		if (selectAllMenuItem == null) {
			selectAllMenuItem = new JMenuItem();
			selectAllMenuItem.setText("Select all");
			selectAllMenuItem.setAccelerator(KeyStroke.getKeyStroke('A', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			selectAllMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.selectAllItems();
				}
			});
		}
		return selectAllMenuItem;
	}

	private JMenu getWorkflowsMenu() {
		if (workflowsMenu == null) {
			workflowsMenu = new JMenu();
			workflowsMenu.setText("Workflow");
			workflowsMenu.add(getOpenWorkflowMenuItem());
			workflowsMenu.add(getRecentWorkflowMenu());
			workflowsMenu.add(getOpenRepositoryWorkflowMenu());
			workflowsMenu.addSeparator();
			workflowsMenu.add(getSaveWorkflowMenuItem());
			if (application.isStandalone) {
				workflowsMenu.setEnabled(false);
			}
		}
		return workflowsMenu;
	}

	private JMenuItem getTaskListMenuItem() {
		if (taskListMenuItem == null) {
			taskListMenuItem = new JMenuItem();
			taskListMenuItem.setText("Tasks...");
			taskListMenuItem.setAccelerator(KeyStroke.getKeyStroke('T', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			taskListMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.flipTaskListVisibility(false);
				}
			});
		}
		return taskListMenuItem;
	}

	private JMenuItem getSaveWorkflowMenuItem() {
		if (saveWorkflowMenuItem == null) {
			saveWorkflowMenuItem = new JMenuItem();
			saveWorkflowMenuItem.setText("Save starting from selected...");
			saveWorkflowMenuItem.setEnabled(false);
			saveWorkflowMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					File workflow = application.saveWorkflow();
					addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));
				}
			});
		}
		return saveWorkflowMenuItem;
	}

	private JMenu getRecentWorkflowMenu() {
		if (recentWorkflowMenu == null) {
			recentWorkflowMenu = new JMenu();
			recentWorkflowMenu.setText("Run recent");

			List<File> workflows = application.getWorkflows();

			for (File workflow : workflows) {
				addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));
			}
		}
		return recentWorkflowMenu;
	}

	private JMenuItem createRunWorkflowMenuItem(String name, final URL workflowScript) {
		JMenuItem runWorkflowMenuItem = new JMenuItem(name);
		runWorkflowMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				application.runWorkflow(workflowScript);
			}
		});
		return runWorkflowMenuItem;
	}


	private JMenu getViewMenu() {
		if (viewMenu == null) {
			viewMenu = new JMenu();
			viewMenu.setText("View");
			viewMenu.add(getRestoreViewMenuItem());
			viewMenu.addSeparator();
			viewMenu.add(getVisualiseMenutItem());
			viewMenu.add(getVisualisationwMenu());
			viewMenu.addSeparator();
			viewMenu.add(getFontSize());
			viewMenu.addSeparator();
			viewMenu.add(getTaskListMenuItem());

		}
		return viewMenu;
	}

	private JMenu getVisualisationwMenu() {
		if (visualisationMenu == null) {
			visualisationMenu = new JMenu();
			visualisationMenu.setText("Visualisation");
			visualisationMenu.add(getMaximiseVisualisationMenuItem());
			visualisationMenu.add(getDetachMenuItem());

		}
		return visualisationMenu;
	}

	private JMenuItem getDetachMenuItem() {
		if (detachMenuItem == null) {
			detachMenuItem = new JMenuItem();
			detachMenuItem.setText("Detach");
			detachMenuItem.setAccelerator(KeyStroke.getKeyStroke('N', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			detachMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					VisualisationToolBar visToolBar = application.getVisualisationFrameManager().getVisualisationToolBar();

					visToolBar.detach();
				}

			});
		}

		return detachMenuItem;
	}

	private JMenuItem getMaximiseVisualisationMenuItem() {
		if (maximiseVisualisationMenuItem == null) {
			maximiseVisualisationMenuItem = new JMenuItem();
			maximiseVisualisationMenuItem.setText("Maximise/Restore");
			maximiseVisualisationMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F4, 0));
			maximiseVisualisationMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					VisualisationToolBar visToolBar = application.getVisualisationFrameManager().getVisualisationToolBar();

					visToolBar.maximiseOrRestoreVisualisation();
				}

			});
		}

		return maximiseVisualisationMenuItem;
	}

	private JMenuItem getVisualiseMenutItem() {
		if (visualiseMenuItem == null) {
			visualiseMenuItem = new JMenuItem();
			visualiseMenuItem.setText("Visualise selected");
			visualiseMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0));
			visualiseMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					application.visualiseWithBestMethod(FrameType.MAIN);
				}

			});
		}
		return visualiseMenuItem;
	}

	private JMenuItem getRestoreViewMenuItem() {
		if (restoreViewMenuItem == null) {
			restoreViewMenuItem = new JMenuItem();
			restoreViewMenuItem.setText("Restore default");
			restoreViewMenuItem.setIcon(VisualConstants.DEFAULT_VIEW_MENUICON);
			restoreViewMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.restoreDefaultView();
				}
			});
		}
		return restoreViewMenuItem;
	}

	/**
	 * This method initializes helpInfoMenu. Name was chosen because getHelpMenu() causes conflict.
	 */
	private JMenu getHelpInfoMenu() {
		if (helpInfoMenu == null) {
			helpInfoMenu = new JMenu();
			helpInfoMenu.setText("Help");
			helpInfoMenu.setMnemonic('H');
			helpInfoMenu.add(getStartedMenuItem());
			helpInfoMenu.add(getContentMenuItem());
			helpInfoMenu.add(getAboutMenuItem());
		}
		return helpInfoMenu;
	}


	private JMenuItem getContentMenuItem() {
		if (contentMenuItem == null) {
			contentMenuItem = new JMenuItem();
			contentMenuItem.setText("User manual");
			contentMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0));
			contentMenuItem.setIcon(VisualConstants.HELP_MENUICON);
			contentMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewHelp("chipster-manual/index.html");
				}
			});
		}
		return contentMenuItem;
	}

	private JMenuItem getStartedMenuItem() {
		if (startedMenuItem == null) {
			startedMenuItem = new JMenuItem();
			startedMenuItem.setText("Getting started");
			startedMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewHelp("chipster-manual/basic-functionality.html");
				}
			});
		}
		return startedMenuItem;
	}

	private JMenuItem getAboutMenuItem() {
		if (aboutMenuItem == null) {
			aboutMenuItem = new JMenuItem();
			aboutMenuItem.setText("About");
			aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewHelp("chipster-manual/about.html");
				}
			});
		}
		return aboutMenuItem;
	}

	private JMenu getFontSize() {
		if (fontSizeMenu == null) {
			fontSizeMenu = new JMenu();
			fontSizeMenu.setText("Text size");

			// FIXME L&F Theme is lost when the font is changed
			// fontSizeMenu.setEnabled(false);

			JMenuItem norm = new JMenuItem("Normal");
			JMenuItem inc = new JMenuItem("Increase");
			JMenuItem dec = new JMenuItem("Decrease");

			inc.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_PLUS, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			dec.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_MINUS, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			norm.setAccelerator(KeyStroke.getKeyStroke('0', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));

			fontSizeMenu.add(inc);
			fontSizeMenu.add(dec);
			fontSizeMenu.addSeparator();
			fontSizeMenu.add(norm);

			norm.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.setFontSize(VisualConstants.DEFAULT_FONT_SIZE);
				}
			});

			dec.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.setFontSize(application.getFontSize() - 1);
				}
			});

			inc.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.setFontSize(application.getFontSize() + 1);
				}
			});
		}
		return fontSizeMenu;
	}

	private JMenuItem getClearSessionMenuItem() {
		if (clearSessionMenuItem == null) {
			clearSessionMenuItem = new JMenuItem();
			clearSessionMenuItem.setText("New session");
			clearSessionMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					application.clearSession();
				}

			});
		}
		return clearSessionMenuItem;
	}

	private JMenuItem getLoadSnapshotMenuItem() {
		if (loadSnapshotMenuItem == null) {
			loadSnapshotMenuItem = new JMenuItem();
			loadSnapshotMenuItem.setText("Open session...");
			loadSnapshotMenuItem.setAccelerator(KeyStroke.getKeyStroke('O', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			loadSnapshotMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						application.loadSession();

					} catch (Exception ioe) {
						application.reportException(ioe);
					}
				}
			});
		}
		return loadSnapshotMenuItem;
	}

	private JMenuItem getSaveSnapshotMenuItem() {
		if (saveSnapshotMenuItem == null) {
			saveSnapshotMenuItem = new JMenuItem();
			saveSnapshotMenuItem.setText("Save session...");
			saveSnapshotMenuItem.setAccelerator(KeyStroke.getKeyStroke('S', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			saveSnapshotMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.saveSession();
				}
			});
		}
		return saveSnapshotMenuItem;
	}

	private JMenuItem getAddDirMenuItem() {
		if (addDirMenuItem == null) {
			addDirMenuItem = new JMenuItem();
			addDirMenuItem.setText("Import folder...");
			addDirMenuItem.setAccelerator(KeyStroke.getKeyStroke('I', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask() | InputEvent.SHIFT_DOWN_MASK, false));
			addDirMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.openDirectoryImportDialog();
				}
			});
		}
		return addDirMenuItem;
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof DatasetChoiceEvent || evt instanceof VisualisationMethodChangedEvent) {

			updateMenuStatus();
		}

	}
}
