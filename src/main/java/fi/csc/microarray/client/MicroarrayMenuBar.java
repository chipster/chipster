package fi.csc.microarray.client;

import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.SwingClientApplication.SessionSavingMethod;
import fi.csc.microarray.client.dialog.ClipboardImportDialog;
import fi.csc.microarray.client.dialog.FeedbackDialog;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationToolBar;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;
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
	private JMenuItem importFromURLToServerMenuItem = null;
	private JMenuItem importFromClipboardMenuItem = null;
	private JMenuItem openWorkflowsMenuItem = null;
	private JMenuItem openWorkflowsForEachMenuItem = null;
	private JMenu recentWorkflowMenu;
	private JMenu recentWorkflowForEachMenu;
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
	private JMenuItem sendFeedbackMenuItem;
	private JMenuItem saveWorkflowMenuItem;
	private JMenuItem helpWorkflowMenuItem;
	private JMenuItem archiveSessionMenuItem;
	private JMenuItem saveSessionMenuItem;
	private JMenuItem taskListMenuItem;
	private JMenuItem clearSessionMenuItem;
	private JMenu joinSessionMenu;
	private JMenuItem selectAllMenuItem;
	private JMenuItem historyMenuItem;
	private JMenuItem maximiseVisualisationMenuItem;
	private JMenuItem visualiseMenuItem;
	private JMenuItem closeVisualisationMenuItem;
	private JMenuItem detachMenuItem;
	private JMenu visualisationMenu;
	private JMenu openRepoWorkflowsMenu;

	private boolean hasRepoWorkflows;

	private JMenuItem manageSessionsMenuItem;

	private JMenuItem openExampleSessionMenuItem;

	public MicroarrayMenuBar(SwingClientApplication application) {
		this.application = application;
		add(getFileMenu());
		add(getEditMenu());
		add(getViewMenu());
		add(getWorkflowsMenu());
		add(getHelpInfoMenu());
		application.addClientEventListener(this);
	}

	public void updateMenuStatus() {

		logger.debug("updating menubar when selected is " + application.getSelectionManager().getSelectedItem());
		DataSelectionManager selectionManager = application.getSelectionManager();

		DataBean selectedDataBean = selectionManager.getSelectedDataBean();
		boolean somethingSelected = selectionManager.getSelectedItem() != null;
		boolean multipleDatasSelected = selectionManager.getSelectedDataBeans().size() > 1;
		
		boolean workflowCompatibleDataSelected = false;

		if (selectedDataBean != null) {
			workflowCompatibleDataSelected = Session.getSession().getPrimaryModule().isWorkflowCompatible(selectedDataBean);
		}
		
		renameMenuItem.setEnabled(selectionManager.getSelectedDataBeans().size() == 1);

		historyMenuItem.setEnabled(selectedDataBean != null && application.getSelectionManager().getSelectedDataBeans().size() == 1);

		visualiseMenuItem.setEnabled(selectedDataBean != null);

		VisualisationMethod method = application.getVisualisationFrameManager().getFrame(FrameType.MAIN).getMethod();
		visualisationMenu.setEnabled(method != null);
		
		closeVisualisationMenuItem.setEnabled(!VisualisationMethod.isDefault(method));

		openWorkflowsMenuItem.setEnabled(workflowCompatibleDataSelected);
		recentWorkflowMenu.setEnabled(workflowCompatibleDataSelected);
		openRepoWorkflowsMenu.setEnabled(workflowCompatibleDataSelected && hasRepoWorkflows);
		openWorkflowsForEachMenuItem.setEnabled(workflowCompatibleDataSelected && multipleDatasSelected);
		recentWorkflowForEachMenu.setEnabled(workflowCompatibleDataSelected && multipleDatasSelected);
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
			fileMenu.add(getLoadLocalSessionMenuItem(true));
			fileMenu.add(getSaveLocalSessionMenuItem());
			fileMenu.addSeparator();
			fileMenu.add(getOpenExampleSessionMenuItem());
			if (application.getSessionManager().areCloudSessionsEnabled()) {
				fileMenu.addSeparator();
				fileMenu.add(getLoadSessionMenuItem(true));
				fileMenu.add(getSaveSessionMenuItem());
				fileMenu.add(getManageSessionsMenuItem());
			}
			fileMenu.addSeparator();			
			fileMenu.add(getMergeSessionMenu());
			fileMenu.add(getClearSessionMenuItem());
			fileMenu.addSeparator();
			fileMenu.add(getQuitMenuItem());
		}
		return fileMenu;
	}

	private JMenuItem getMergeSessionMenu() {
		if (joinSessionMenu == null) {
			joinSessionMenu = new JMenu();
			joinSessionMenu.setText("Merge session");
			joinSessionMenu.add(getLoadLocalSessionMenuItem(false));
			if (application.getSessionManager().areCloudSessionsEnabled()) {
				joinSessionMenu.add(getLoadSessionMenuItem(false));
			}
		}
		return joinSessionMenu;
	}

	private JMenu getImportMenu() {
		if (importMenu == null) {
			importMenu = new JMenu();
			importMenu.setText("Import from");
			
			if (!application.isStandalone()) {
				Module primaryModule = Session.getSession().getPrimaryModule();
				primaryModule.addImportMenuItems(importMenu);
			}
				
			importMenu.add(getImportFromURLMenuItem());
			importMenu.add(getImportFromURLToServerMenuItem());
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
			importFromURLMenuItem.setText("URL to client...");
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

	private JMenuItem getImportFromURLToServerMenuItem() {
		if (importFromURLToServerMenuItem == null) {
			importFromURLToServerMenuItem = new JMenuItem();
			importFromURLToServerMenuItem.setText("URL directly to server...");
			importFromURLToServerMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					BasicModule.importFromUrlToServer();
				}
			});
		}
		return importFromURLToServerMenuItem;
	}

	private JMenuItem getHelpWorkflowMenuItem() {
		if (helpWorkflowMenuItem == null) {
			helpWorkflowMenuItem = new JMenuItem();
			helpWorkflowMenuItem.setText("More information...");
			helpWorkflowMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						application.viewHelp("workflows.html#ready");
						
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
						File workflow = application.openWorkflow(false);
						if (workflow != null) {
							addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));
						}
					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return openWorkflowsMenuItem;
	}

	private JMenuItem getOpenWorkflowForEachMenuItem() {
		if (openWorkflowsForEachMenuItem == null) {
			openWorkflowsForEachMenuItem = new JMenuItem();
			openWorkflowsForEachMenuItem.setText("Run for each...");
			openWorkflowsForEachMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					try {
						File workflow = application.openWorkflow(true);
						addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));

					} catch (Exception me) {
						application.reportException(me);
					}
				}
			});
		}
		return openWorkflowsForEachMenuItem;
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
			((SwingClientApplication) application).runWorkflow(url, false);
			addRecentWorkflow(name, url);
		}
	}

	public void addRecentWorkflow(String name, URL url) {
		if (url != null) {
			// Check if this exists already
			for (int i = 0; i < recentWorkflowMenu.getItemCount(); i++) {
				JMenuItem menuItem = recentWorkflowMenu.getItem(i);
				if (menuItem.getText().equals(name)) {
					recentWorkflowMenu.remove(menuItem);
				}

			}
			
			// Update both recent lists (normal and for-each)
			recentWorkflowMenu.add(createRunWorkflowMenuItem(name, url, false));
			recentWorkflowForEachMenu.add(createRunWorkflowMenuItem(name, url, true));
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
			historyMenuItem.setIcon(VisualConstants.getIcon(VisualConstants.GENERATE_HISTORY_ICON));
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
	
					application.showRenameView();
				}
			});
		}
		return renameMenuItem;
	}

	private JMenuItem getDeleteMenuItem() {
		if (deleteMenuItem == null) {
			deleteMenuItem = new JMenuItem();
			deleteMenuItem.setText("Delete selected item");
			deleteMenuItem.setIcon(VisualConstants.getIcon(VisualConstants.DELETE_MENUICON));
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
			workflowsMenu.add(getOpenWorkflowForEachMenuItem());
			workflowsMenu.add(getRecentWorkflowForEachMenu());
			workflowsMenu.addSeparator();
			workflowsMenu.add(getSaveWorkflowMenuItem());

			if (application.isStandalone()) {
				workflowsMenu.setEnabled(false);
				
			} else {
				
				// Populate both recent workflow lists
				List<File> workflows = application.getWorkflows();
				for (File workflow : workflows) {
					addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));
				}
			}
		}
		return workflowsMenu;
	}

	private JMenuItem getTaskListMenuItem() {
		if (taskListMenuItem == null) {
			taskListMenuItem = new JMenuItem();
			taskListMenuItem.setText("View jobs...");
			taskListMenuItem.setAccelerator(KeyStroke.getKeyStroke('T', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			taskListMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewTasks();
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
					if (workflow != null) {
						addRecentWorkflow(workflow.getName(), Files.toUrl(workflow));
					}
				}
			});
		}
		return saveWorkflowMenuItem;
	}

	private JMenu getRecentWorkflowMenu() {
		if (recentWorkflowMenu == null) {
			recentWorkflowMenu = new JMenu();
			recentWorkflowMenu.setText("Run recent");
		}
		return recentWorkflowMenu;
	}

	private JMenu getRecentWorkflowForEachMenu() {
		if (recentWorkflowForEachMenu == null) {
			recentWorkflowForEachMenu = new JMenu();
			recentWorkflowForEachMenu.setText("Run recent for each");
		}
		return recentWorkflowForEachMenu;
	}

	private JMenuItem createRunWorkflowMenuItem(String name, final URL workflowScript, final boolean runForEach) {
		JMenuItem runWorkflowMenuItem = new JMenuItem(name);
		runWorkflowMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				application.runWorkflow(workflowScript, runForEach);
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
			viewMenu.add(getCloseVisualisationMenutItem());
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
	
	private JMenuItem getCloseVisualisationMenutItem() {
		if (closeVisualisationMenuItem == null) {
			closeVisualisationMenuItem = new JMenuItem();
			closeVisualisationMenuItem.setText("Close visualisation");
			closeVisualisationMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0));
			closeVisualisationMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					application.setVisualisationMethodToDefault();
				}

			});
		}
		return closeVisualisationMenuItem;
	}

	private JMenuItem getRestoreViewMenuItem() {
		if (restoreViewMenuItem == null) {
			restoreViewMenuItem = new JMenuItem();
			restoreViewMenuItem.setText("Restore default");
			restoreViewMenuItem.setIcon(VisualConstants.getIcon(VisualConstants.DEFAULT_VIEW_MENUICON));
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
			if (DirectoryLayout.getInstance().getConfiguration().getBoolean("client", "enable-contact-support")) {
				helpInfoMenu.add(getSendFeedbackMenuItem());
			}

			helpInfoMenu.add(getAboutMenuItem());
		}
		return helpInfoMenu;
	}


	private JMenuItem getContentMenuItem() {
		if (contentMenuItem == null) {
			contentMenuItem = new JMenuItem();
			contentMenuItem.setText("User manual");
			contentMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0));
			contentMenuItem.setIcon(VisualConstants.getIcon(VisualConstants.HELP_MENUICON));
			contentMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewHelp(Session.getSession().getPrimaryModule().getManualHome());
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
					application.viewHelp(Session.getSession().getPrimaryModule().getManualHome() + "/basic-functionality.html");
				}
			});
		}
		return startedMenuItem;
	}

	private JMenuItem getSendFeedbackMenuItem() {
		if (sendFeedbackMenuItem == null) {
			sendFeedbackMenuItem = new JMenuItem();
			sendFeedbackMenuItem.setText("Contact Support...");
			sendFeedbackMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					FeedbackDialog feedback = new FeedbackDialog(application, "", false);
                    feedback.showDialog();
				}
			});
		}
		return sendFeedbackMenuItem;
	}
	
	
	private JMenuItem getAboutMenuItem() {
		if (aboutMenuItem == null) {
			aboutMenuItem = new JMenuItem();
			aboutMenuItem.setText("About");
			aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.viewHelp(Session.getSession().getPrimaryModule().getManualHome() + "/about.html");
				}
			});
		}
		return aboutMenuItem;
	}

	private JMenu getFontSize() {
		if (fontSizeMenu == null) {
			fontSizeMenu = new JMenu();
			fontSizeMenu.setText("Text size");

			// TODO fix: L&F Theme is lost when the font is changed
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
	
	private JMenuItem getOpenExampleSessionMenuItem() {
		if (openExampleSessionMenuItem == null) {
			openExampleSessionMenuItem = new JMenuItem();
			openExampleSessionMenuItem.setText("Open example session");
			openExampleSessionMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					application.loadSession(true, true, true);
				}

			});
		}
		return openExampleSessionMenuItem;
	}

	private JMenuItem getClearSessionMenuItem() {
		if (clearSessionMenuItem == null) {
			clearSessionMenuItem = new JMenuItem();
			clearSessionMenuItem.setText("New session");
			clearSessionMenuItem.addActionListener(new ActionListener() {

				public void actionPerformed(ActionEvent e) {
					try {
						application.clearSession();
					} catch (MalformedURLException | FileBrokerException | AuthCancelledException e1) {
						application.reportException(e1);
					}
				}

			});
		}
		return clearSessionMenuItem;
	}

	private JMenuItem getLoadSessionMenuItem(final boolean clear) {
		JMenuItem loadSessionMenuItem = new JMenuItem();
		if (clear) {			
			loadSessionMenuItem.setText("Open cloud session... (BETA)");
		} else {
			loadSessionMenuItem.setText("Cloud session... (BETA)");
		}
		loadSessionMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				try {
					application.loadSession(true, false, clear);

				} catch (Exception ioe) {
					application.reportException(ioe);
				}
			}
		});
		return loadSessionMenuItem;
	}

	private JMenuItem getLoadLocalSessionMenuItem(final boolean clear) {
		JMenuItem loadLocalSessionMenuItem = new JMenuItem();
		loadLocalSessionMenuItem.setText("Open local session...");
		if (clear) {
			loadLocalSessionMenuItem.setText("Open local session...");
			loadLocalSessionMenuItem.setAccelerator(KeyStroke.getKeyStroke('O', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
		} else {
			loadLocalSessionMenuItem.setText("Local session...");
		}
		loadLocalSessionMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				try {
					application.loadSession(false, false, clear);
				} catch (Exception ioe) {
					application.reportException(ioe);
				}
			}
		});
		return loadLocalSessionMenuItem;
	}

	private JMenuItem getSaveSessionMenuItem() {
		if (saveSessionMenuItem == null) {
			saveSessionMenuItem = new JMenuItem();
			saveSessionMenuItem.setText("Save cloud session... (BETA)");
			saveSessionMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.saveSession(SessionSavingMethod.UPLOAD_DATA_TO_SERVER);
				}
			});
		}
		return saveSessionMenuItem;
	}

	private JMenuItem getManageSessionsMenuItem() {
		if (manageSessionsMenuItem == null) {
			manageSessionsMenuItem = new JMenuItem();
			manageSessionsMenuItem.setText("Manage cloud sessions... (BETA)");
			manageSessionsMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.manageRemoteSessions();
				}
			});
		}
		return manageSessionsMenuItem;
	}


	private JMenuItem getSaveLocalSessionMenuItem() {
		if (archiveSessionMenuItem == null) {
			archiveSessionMenuItem = new JMenuItem();
			archiveSessionMenuItem.setText("Save local session...");
			archiveSessionMenuItem.setAccelerator(KeyStroke.getKeyStroke('S', Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(), false));
			archiveSessionMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					application.saveSession(SessionSavingMethod.INCLUDE_DATA_INTO_ZIP);
				}
			});
		}
		return archiveSessionMenuItem;
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
