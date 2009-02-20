package fi.csc.microarray.client.dataimport;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileSystemView;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.table.InformationDialog;
import fi.csc.microarray.client.dialog.ChipsterDialog;
import fi.csc.microarray.client.dialog.DialogInfo;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.dialog.DialogInfo.Type;
import fi.csc.microarray.config.MicroarrayConfiguration;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.util.IOUtils;

/**
 * Util class for the import dataset choosers (JFileChooser, URL import,
 * clipboard paste). Contains methods to help implementation of folder selection
 * and launching actionChooser or direct import.
 * 
 * @author Petri Klemelä
 */
public class ImportUtils {

	private static final Logger logger = Logger.getLogger(ImportUtils.class);
	private static final String DEFAULT_FOLDER_NAME = "My experiment";
	private static ClientApplication application = Session.getSession().getApplication();

	private static boolean zipDialogShown = false;

	/**
	 * <strong>You must always use this!</strong> This is a convenience method,
	 * see the actual for details.
	 * 
	 * @see #getFixedFileChooser(File)
	 */
	public static JFileChooser getFixedFileChooser() {
		return getFixedFileChooser(null);
	}

	/**
	 * <p>
	 * Returns a file chooser that should not trip the bug with Java and Windows
	 * XP zip/rar folder feature. Also checks and reports with a modal dialog if
	 * there are zips on desktop, as the does not seem to work in all
	 * environments. For more information about the bug see e.g. #5050516 at
	 * bugs.sun.com.
	 * </p>
	 * 
	 * <p>
	 * <strong>You must always use this!</strong>
	 * </p>
	 * 
	 * @return a safe file chooser
	 */
	public static JFileChooser getFixedFileChooser(File file) {

		// first check for the bug
		String osName = System.getProperty("os.name");
		boolean hasZipFolderSupport = osName.contains("Windows XP") || osName.contains("Windows NT (unknown)") || osName.contains("Windows Vista");

		if (hasZipFolderSupport && !zipDialogShown) {
			boolean zipsOnDesktop = false;

			// check desktop folder
			File[] roots = FileSystemView.getFileSystemView().getRoots();
			for (File root : roots) {

				// list zips
				String[] zips = root.list(new FilenameFilter() {
					public boolean accept(File dir, String name) {
						String lcName = name.toLowerCase();
						return lcName.endsWith(".zip");
					}
				});

				if (zips.length > 0) {
					zipsOnDesktop = true;
					break;
				}
			}

			if (zipsOnDesktop) {
				DialogInfo info = new DialogInfo(Severity.INFO, "ZIP files detected", 
						"There seems to be ZIP files on your desktop. "
						+ "This can slow down selecting files in some Windows versions. " 
						+ "If you experience this problem, please move the ZIP files to a subfolder.", null, Type.OK_MESSAGE); 
				ChipsterDialog.showDialog(null, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true, null);
				zipDialogShown = true;
			}
		}

		JFileChooser fileChooser = file != null ? new JFileChooser(file) : new JFileChooser();
		fileChooser.putClientProperty("FileChooser.useShellFolder", Boolean.FALSE);

		return fileChooser;
	}

	/**
	 * @return Set of Strings containing names of dataset folders from the
	 *         client
	 */
	public static Set<String> getFolderNames(boolean includeDefaultFolder) {
		Set<String> folderNameList = new HashSet<String>();
		DataFolder root = application.getDataManager().getRootFolder();

		// TODO Is it necessary to put file into root Folder?
		// folderNameList.add(root.getName());
		for (DataItem item : root.getChildren()) {
			if (item instanceof DataFolder) {
				folderNameList.add(item.getName());
			}
		}

		if (includeDefaultFolder) {
			folderNameList.add(DEFAULT_FOLDER_NAME);
		}

		return folderNameList;
	}

	public static ImportUtils.URLFileLoader getURLFileLoader() {
		return new ImportUtils.URLFileLoader();
	}

	public static class URLFileLoader {

		public File loadFileFromURL(URL url, File outputFile, String importFolder, boolean skipActionChooser) {
			logger.debug("Method loadFileFromURL started");
			InformationDialog info = new InformationDialog("Loading file", "Loading file from the specified URL", null);
			logger.debug("Next the download process will start");
			new FileLoaderProcess(outputFile, url, importFolder, info, skipActionChooser).runProcess();
			logger.debug("Download process started");
			return outputFile;
		}
	}

	public static class FileLoaderProcess extends RunnableImportProcess {

		private File outputFile;
		private URL url;
		private InformationDialog info;
		private String importFolder;
		private boolean skipActionChooser;

		public FileLoaderProcess(File outputFile, URL url, String importFolder, InformationDialog info, boolean skipActionChooser) {
			super(info);
			this.outputFile = outputFile;
			this.url = url;
			this.info = info;
			this.importFolder = importFolder;
			this.skipActionChooser = skipActionChooser;
		}

		public void taskToDo() {

			HttpURLConnection connection = null;
			try {
				connection = (HttpURLConnection) url.openConnection();
				info.setMinimunValue(0);
				info.setMaximumValue(connection.getContentLength());

				BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
				InputStream in = connection.getInputStream();
				byte[] buffer = new byte[1024];
				int numRead;
				int numWritten = 0;
				while ((numRead = in.read(buffer)) != -1) {
					out.write(buffer, 0, numRead);
					numWritten += numRead;
					info.setValue(numWritten);
				}

				// This is important, else small files are lost
				out.close();

				if (numWritten > 0) {
					SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							ImportUtils.executeImport(new ImportSession(ImportSession.Source.URL, new File[] { outputFile }, importFolder, skipActionChooser));
						}
					});
				} else {
					JOptionPane.showMessageDialog(((SwingClientApplication)application).getMainFrame(), "Length of the loaded file is zero, import aborted", "File size too small", JOptionPane.ERROR_MESSAGE);
				}
				
			} catch (IOException e) {
				if (e instanceof FileNotFoundException) {
					JOptionPane.showMessageDialog(((SwingClientApplication) application).getMainFrame(), "File from the typed url can't be found", "File not found", JOptionPane.ERROR_MESSAGE);
				} else if (e instanceof UnknownHostException) {
					JOptionPane.showMessageDialog(((SwingClientApplication) application).getMainFrame(), "Host from the typed url can't be found", "Host not found", JOptionPane.ERROR_MESSAGE);
				} else {
					application.reportException(e);
				}
				
			} catch (IllegalArgumentException e) {
				JOptionPane.showMessageDialog(((SwingClientApplication) application).getMainFrame(), "Malformed URL, correct URL form is http://www.host.com/file.ext", "Malformed URL", JOptionPane.ERROR_MESSAGE);
				
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
			
		}
	}

	/**
	 * Creates a new temp file, which is going to be removed when the JVM is
	 * closed. Name is used as a start of the filename, put some extra
	 * characters will be appended to it.
	 * 
	 * @param name
	 * @return
	 * @throws IOException
	 */
	public static File createTempFile(String name, String ext) throws IOException {
		File file = File.createTempFile(name, ext);
		file.deleteOnExit();
		return file;
	}

	/**
	 * Iterates through the files and returns true if any of the files is not
	 * accepted with the isFileSupported() method.
	 * 
	 * @param files
	 * @return false if all files are supported type
	 */
	public static boolean containsUnsupportedTypes(File[] files) {
		boolean isUnsupported = false;

		for (File file : files) {
			if (!ImportUtils.isFileSupported(file)) {
				isUnsupported = true;
			}
		}

		return isUnsupported;
	}

	/**
	 * @param file
	 * @return true if Microarray Files -filefilter accepts the given file
	 */
	public static boolean isFileSupported(File file) {
		return application.getDataManager().guessContentType(file).isSupported();
	}

	/**
	 * @param fileName
	 * @return fileName without file type extension
	 */
	// TODO Remove this or implements it in some other way. Splitting
	// using dot is not the best way to catch the filename extension
	public static String convertToDatasetName(String fileName) {
		// for now, do nothing
		return fileName;
		/*
		 * if(fileName.contains(".")){ return fileName.substring(0,
		 * fileName.indexOf('.')); } else { return fileName; }
		 */
	}

	// TODO Remove this or implements it in some other way. Splitting
	// using dot is not the best way to catch the filename extension
	public static String getExtension(String fileName) {
		if (fileName.contains(".")) {
			return fileName.substring(fileName.indexOf('.'), fileName.length());
		} else {
			return "";
		}
	}

	public static String URLToFilename(URL url) {
		String fileName;
		// Removes the folder (for example /path/to/file.ext -> file.ext)
		fileName = url.getFile().substring(url.getFile().lastIndexOf("/") + 1, url.getFile().length());
		if (ImportUtils.convertToDatasetName(fileName).length() < 3) {
			fileName = "url_" + fileName;
		}

		return fileName;
	}

	/**
	 * Imports given files if they all are supported type or launches
	 * ActionChooserScreen if not.
	 * 
	 * @param files
	 *            as File objects
	 * @param skipActionChooser
	 * @param importFolder
	 *            Folder to where imported data is put
	 */
	public static void executeImport(ImportSession importSession) {

		List<File> files = importSession.getInputFiles();

		if (importSession.isSkipActionChooser() && !ImportUtils.containsUnsupportedTypes(files.toArray(new File[files.size()]))) {
			// skip requested and all of the files are supported => import directly and don't show action chooser			
			application.importGroup(importSession.getImportItems(), importSession.getDestinationFolder());
			
		} else {
			// skip not requested => show ActionChooser
			ActionChooserScreen actionChooser = new ActionChooserScreen(importSession);
			actionChooser.getFrame().setVisible(true);

		}
	}
	
	public static void main(String[] args) throws IOException, OldConfigurationFormatException, MicroarrayException {
		MicroarrayConfiguration.loadConfiguration();
		
		DialogInfo info = new DialogInfo(Severity.INFO, "ZIP files detected", 
				"There seems to be ZIP files on your desktop. "
				+ "This can slow down selecting files in some Windows versions. " 
				+ "If you experience this problem, please move the ZIP files to a subfolder.", null, Type.OK_MESSAGE); 
		ChipsterDialog.showDialog(null, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true, null);

	}
}
