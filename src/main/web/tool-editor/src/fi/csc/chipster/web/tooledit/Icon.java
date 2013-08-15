package fi.csc.chipster.web.tooledit;

import java.io.File;

import com.vaadin.server.FileResource;
import com.vaadin.server.Resource;
import com.vaadin.server.VaadinService;

/**
 * Class for getting icons
 * @author Gintare Pacauskaite
 *
 */
public class Icon {
	
	private static final String[] pathEnd = new String[] {"web", "tool-editor", "WebContent"};
	private static final String separator = File.separator;

	public static String getUpButtonIconPath() {
		String[] deletePath = new String[] {"resources", "arrow_up.png"};
		return removePathEnd() + getGeneratedPath(deletePath);
	}
	
	public static String getDownButtonIconPath() {
		String[] deletePath = new String[] {"resources", "arrow_down.png"};
		return removePathEnd() + getGeneratedPath(deletePath);
	}
	
	public static String getDeleteButtonIconPath() {
		String[] deletePath = new String[] {"resources", "close.png"};
		return removePathEnd() + getGeneratedPath(deletePath);
	}
	
	public static String getAddButtonIconPath() {
		String[] deletePath = new String[] {"resources", "edit-add-2.png"};
		return removePathEnd() + getGeneratedPath(deletePath);
	}
	
	private static String removePathEnd() {
		return VaadinService.getCurrent().getBaseDirectory().getAbsolutePath().replace(getPathEnd(), "");
	}
	
	private static String getPathEnd() {
		return getGeneratedPath(pathEnd);
	}
	
	private static String getGeneratedPath(String[] array) {
		StringBuilder string = new StringBuilder();
		for(int i = 0; i < array.length; i++) {
			string.append(separator);
			string.append(array[i]);
		}
		return string.toString();
	}
	
	public static Resource getResource(String path) {
		return new FileResource(new File(path));
	}

}
