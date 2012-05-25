package fi.csc.microarray.module.basic;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collections;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.filechooser.FileFilter;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.ExternalBrowserViewer;
import fi.csc.microarray.client.visualisation.methods.HtmlViewer;
import fi.csc.microarray.client.visualisation.methods.ImageViewer;
import fi.csc.microarray.client.visualisation.methods.PDFViewer;
import fi.csc.microarray.client.visualisation.methods.Spreadsheet;
import fi.csc.microarray.client.visualisation.methods.TextViewer;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.TypeTag;
import fi.csc.microarray.databeans.features.RestrictModifier;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.databeans.features.bio.PhenodataProvider;
import fi.csc.microarray.databeans.features.stat.LogModifier;
import fi.csc.microarray.databeans.features.stat.NegModifier;
import fi.csc.microarray.databeans.features.table.RowCountProvider;
import fi.csc.microarray.databeans.features.table.TableColumnProvider;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.util.GeneralFileFilter;

public class BasicModule implements Module {

	public static class TypeTags {
		public static final TypeTag TABLE_WITHOUT_COLUMN_NAMES = new TypeTag("table-without-column-names", "first row is the first data row");
		public static final TypeTag TABLE_WITH_COLUMN_NAMES = new TypeTag("table-without-column-names", "first row is the column name row");
		public static final TypeTag TABLE_WITH_TITLE_ROW = new TypeTag("table-with-possible-title-row", "first row is title row");
		public static final TypeTag PHENODATA = new TypeTag("phenodata", "phenodata table");
	}
	
	public static class VisualisationMethods {
		public static VisualisationMethod SPREADSHEET = new VisualisationMethod("Spreadsheet", Spreadsheet.class, VisualConstants.SPREADSHEET_MENUICON, 2, 0.0007);
		public static VisualisationMethod SHOW_IMAGE = new VisualisationMethod("Show image", ImageViewer.class, VisualConstants.IMAGE_MENUICON, 1, 0.015); 
		public static VisualisationMethod WEBVIEW = new VisualisationMethod("View page", HtmlViewer.class, VisualConstants.HTML_MENUICON, 1, 0.008); 
		public static VisualisationMethod PDFVIEW = new VisualisationMethod("View PDF", PDFViewer.class, VisualConstants.IMAGE_MENUICON, 1, 0);
		public static VisualisationMethod VIEW_TEXT = new VisualisationMethod("View text", TextViewer.class, VisualConstants.TEXT_MENUICON, 1, 0);
		public static VisualisationMethod EXTERNAL_BROWSER = new VisualisationMethod("Open in external web browser", ExternalBrowserViewer.class, VisualConstants.EMPTY_MENUICON, -1, -1);
	}
	
	public void plugContentTypes(DataManager manager) {

		manager.plugContentType("text/plain", false, false, "plain text", VisualConstants.ICON_TYPE_TEXT, "txt", "dat", "wee", "seq", "log", "sam", "fastq");
		manager.plugContentType("application/octet-stream", false, true, "binary", VisualConstants.ICON_TYPE_BINARY, "");
		
		manager.plugContentType("text/tab", false, false, "tab separated values", VisualConstants.ICON_TYPE_TABLE, "tsv");
		manager.plugContentType("text/csv", false, false, "comma separated values", VisualConstants.ICON_TYPE_TABLE, "csv");

		manager.plugContentType("image/png", true, true, "PNG image", VisualConstants.ICON_TYPE_IMAGE, "png");
		manager.plugContentType("image/gif", true, true, "GIF image", VisualConstants.ICON_TYPE_IMAGE, "gif");
		manager.plugContentType("image/jpeg", true, true, "JPEG image", VisualConstants.ICON_TYPE_IMAGE, "jpg", "jpeg");
		
		manager.plugContentType("application/pdf", true, true, "PDF document", VisualConstants.ICON_TYPE_IMAGE, "pdf");		

		manager.plugContentType("text/html", true, false, "HTML document", VisualConstants.ICON_TYPE_HTML, "html", "htm");
		manager.plugContentType("application/pdf", true, true, "PDF document", VisualConstants.ICON_TYPE_IMAGE, "pdf");		
	}

	@Override
	public String[] getServerModuleNames() {
		return null;
	}

	@Override
	public String getModuleLongName(String moduleName) {
		return null;
	}

	@Override
	public void addImportMenuItems(JMenu importMenu) {
		// do nothing
	}

	@Override
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks) {
		// do nothing
	}

	@Override
	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/phenodata", new PhenodataProvider()); // FIXME should be in microarray module, but phenodata checks must be fixed first
		manager.plugFeatureFactory("/column", new TableColumnProvider());
		manager.plugFeatureFactory("/rowcount", new RowCountProvider());
	}

	@Override
	public void plugModifiers(DataManager manager) {
		manager.plugModifier("log", new LogModifier());
		manager.plugModifier("neg", new NegModifier());
		manager.plugModifier("restrict", new RestrictModifier());
	}

	@Override
	public void plugTypeTags(DataManager manager) {
		manager.plugTypeTag(TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
	}

	@Override
	public boolean isImportToolSupported() {
		return false;
	}

	@Override
	public boolean isWorkflowCompatible(DataBean data) {
		return true; // we have to assume that all operations are workflow compatible
	}

	@Override
	public VisualisationMethod[] getVisualisationMethods() {
		
		return new VisualisationMethod[] {
				VisualisationMethod.NONE,
				VisualisationMethods.SPREADSHEET,
				VisualisationMethods.SHOW_IMAGE, 
				VisualisationMethods.WEBVIEW, 
				VisualisationMethods.PDFVIEW,
				VisualisationMethods.VIEW_TEXT,
				VisualisationMethods.EXTERNAL_BROWSER
		};
	}

	@Override
	public URL[] getExampleSessionUrls(boolean isStandalone) throws MalformedURLException {
		return null;
	}

	@Override
	public String[][] getRepositoryWorkflows() {
		return new String[0][0];
	}

	@Override
	public FileFilter[] getImportFileFilter() {
		return new FileFilter[] {
				new GeneralFileFilter("Tab or comma separated tables", new String[] {"csv", "tsv"}),
		};
	}

	@Override
	public boolean isMetadata(DataBean data) {
		return false; // we don't use metadata
	}

	@Override
	public void postProcessOutputMetadata(Operation operation, DataBean metadataOutput) throws MicroarrayException, IOException {
		// do nothing, we don't use metadata
	}

	@Override
	public String getShortCategoryName(String categoryName) {
		return shortenCategoryName(categoryName);
	}

	public static String shortenCategoryName(String catName) {
		if (catName.startsWith("Import")) {
			return "file";
		}
		return catName.substring(0, catName.length() > 4 ? 4 : catName.length());
	}

	@Override
	public boolean countOperationResults() {
		return true;
	}

	@Override
	public JPanel getContextLinkPanel(int selectedDataCount) {
		return null;
	}

	@Override
	public boolean notesVisibleAtStartup() {
		return false;
	}

	@Override
	public String getDisplayName() {
		return "Chipster";
	}

	@Override
	public String getManualHome() {
		return "http://chipster.csc.fi/manual/index.html";
	}

	@Override
	public void addSpeadsheetMenuItems(JPopupMenu spreadsheetPopupMenu, VisualisationFrame visualisationFrame) {
		// do nothing
	}

	@Override
	public List<Boolean> flagLinkableColumns(String[] columnNames) {
		return Collections.nCopies(columnNames.length, false);
	}

	@Override
	public IntegratedEntity createLinkableEntity(Table columns, int column) {
		return null;
	}
}
