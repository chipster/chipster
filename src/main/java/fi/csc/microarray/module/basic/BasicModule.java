package fi.csc.microarray.module.basic;

import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.Module;

public class BasicModule implements Module {

	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("image/png", true, true, "PNG image", VisualConstants.ICON_TYPE_IMAGE, "png");
		manager.plugContentType("image/gif", true, true, "GIF image", VisualConstants.ICON_TYPE_IMAGE, "gif");
		manager.plugContentType("image/jpeg", true, true, "JPEG image", VisualConstants.ICON_TYPE_IMAGE, "jpg", "jpeg");
		
		manager.plugContentType("text/html", true, false, "HTML document", VisualConstants.ICON_TYPE_HTML, "html", "htm");

		manager.plugContentType("text/plain", false, false, "plain text", VisualConstants.ICON_TYPE_TEXT, "txt", "dat", "wee");
		manager.plugContentType("application/octet-stream", false, true, "binary", VisualConstants.ICON_TYPE_BINARY, "");
		
		// FIXME should be separated into different module
		manager.plugContentType("chemical/x-fasta", true, false, "FASTA", VisualConstants.ICON_TYPE_TEXT, "fasta", "fa", "fna", "fsa", "mpfa");
		manager.plugContentType("text/wig", true, false, "WIG", VisualConstants.ICON_TYPE_TEXT, "wig");
		manager.plugContentType("text/bed", true, false, "BED", VisualConstants.ICON_TYPE_TEXT, "bed");
		manager.plugContentType("text/bed-reads", true, false, "READS", VisualConstants.ICON_TYPE_TEXT, "reads");
	}

	public void plugFeatures(DataManager manager) {
		// nothing to plug
	}

	public void plugModifiers(DataManager manager) {
		// nothing to plug
	}

}
