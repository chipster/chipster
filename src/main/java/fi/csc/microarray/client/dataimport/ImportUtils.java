package fi.csc.microarray.client.dataimport;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.table.InformationDialog;
import fi.csc.microarray.client.dialog.TaskImportDialog;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.IOUtils;

/**
 * Util class for the import dataset choosers (JFileChooser, URL import,
 * clipboard paste). Contains methods to help implementation of folder selection
 * and launching actionChooser or direct import.
 * 
 * @author Petri Klemelä
 */
public class ImportUtils {

	private static final String PREPROCESS_TOOL_ID = "LocalNGSPreprocess.java";
	private static final Logger logger = Logger.getLogger(ImportUtils.class);
	private static final String DEFAULT_FOLDER_NAME = "";
	private static ClientApplication application = Session.getSession().getApplication();

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
			new FileLoaderImportProcess(outputFile, url, importFolder, info, skipActionChooser).runProcess();
			logger.debug("Download process started");
			return outputFile;
		}
	}

	public static class FileLoaderProcess extends RunnableImportProcess {

		protected File outputFile;
		protected URL url;
		protected InformationDialog info;

		public FileLoaderProcess(File outputFile, URL url, InformationDialog info) {
			super(info);
			this.outputFile = outputFile;
			this.url = url;
			this.info = info;
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
							postProcess();
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
		
		/**
		 * The default post process does nothing. Subclasses can override to further process copied files.
		 */
		protected void postProcess() {
			// do nothing
		}
	}
	
	public static class FileLoaderImportProcess extends FileLoaderProcess {

		protected String importFolder;
		protected boolean skipActionChooser;

		public FileLoaderImportProcess(File outputFile, URL url, String importFolder, InformationDialog info, boolean skipActionChooser) {
			super(outputFile, url, info);
			this.importFolder = importFolder;
			this.skipActionChooser = skipActionChooser;
		}
		
		@Override
		protected void postProcess() {
			ImportUtils.executeImport(new ImportSession(ImportSession.Source.URL, new File[] { outputFile }, importFolder, skipActionChooser));
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
		if (fileName.length() < 3) {
			fileName = "url_" + fileName;
		}

		return fileName;
	}

	/**
	 * Imports given files if they all are supported type and 
	 * skip is requested (by default it is). Otherwise launches
	 * ActionChooserScreen. If module does not support import
	 * tools then files are always imported directly.
	 * @throws MicroarrayException 
	 */
	public static void executeImport(ImportSession importSession) {

		if (!application.isStandalone()) {
			List<File> files = importSession.getInputFiles();
			
			// bam sam bed gtf and vcf go always to preprocess dialog
			boolean allBamSamOrBed = true;
			for (File file : files) {
				if (!file.getName().toLowerCase().endsWith(".bam") &&
					!file.getName().toLowerCase().endsWith(".bai") &&
					!file.getName().toLowerCase().endsWith(".sam") &&
					!file.getName().toLowerCase().endsWith(".bed") &&
					!file.getName().toLowerCase().endsWith(".gtf") &&
					!file.getName().toLowerCase().endsWith(".vcf")) {
					allBamSamOrBed = false;
					break;
				}
			}

			// go to preprocess with .bam .sam or .bed
			if (allBamSamOrBed) {
				openPreprocessDialog(importSession);
			} 
			
			// normal import
			else {

				boolean importToolSupported = Session.getSession().getPrimaryModule().isImportToolSupported();

				// import directly
				if (!importToolSupported || (importSession.isSkipActionChooser() && !ImportUtils.containsUnsupportedTypes(files.toArray(new File[files.size()])))) {
					// skip requested and all of the files are supported => import directly and don't show action chooser			
					application.importGroup(importSession.getImportItems(), importSession.getDestinationFolder());
				} 

				// action chooser or preprocess
				else {
					new ActionChooserScreen(importSession);
				}
			}
		}

		// standalone
		else {

			// import directly
			if (importSession.isSkipActionChooser()) {
				application.importGroup(importSession.getImportItems(), importSession.getDestinationFolder());
			}

			// go to preprocessing
			else {
				openPreprocessDialog(importSession);
			}
		}
	}

	private static void openPreprocessDialog(ImportSession importSession) {
		// input files to input DataBeans
		try {
			List<DataBean> inputBeans = new LinkedList<DataBean>();
			int i = 0;
			for (File inputFile: importSession.getInputFiles()) {
				inputBeans.add(Session.getSession().getDataManager().createDataBean("preprocessInput-" + i, inputFile));
				i++;
			}

			// create operation, open import operation parameter dialog
			ClientApplication application = Session.getSession().getApplication();
			Operation importOperation = new Operation(application.getOperationDefinition(PREPROCESS_TOOL_ID), inputBeans.toArray(new DataBean[] {}));
			new TaskImportDialog(application, "Preprocess NGS data", importSession, importOperation, "Preprocess", "Cancel", "Skip preprocessing", "Please note that preprocessing SAM and BAM files can take several minutes depending on the file size.");

		} catch (Exception me) {
			Session.getSession().getApplication().reportException(me);
		}
	}	
}
