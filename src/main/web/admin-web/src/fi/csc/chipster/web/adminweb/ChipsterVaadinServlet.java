package fi.csc.chipster.web.adminweb;

import java.net.URL;
import java.util.logging.Logger;

import javax.servlet.http.HttpServletRequest;

import com.vaadin.server.VaadinServlet;

/**
 * Original implementation of this method doesn't allow loading of ThemeResources. ThemeResources
 * have path "!/WebContent/VAADIN/", while only "!/VAADIN/" is allowed in the original implementation.
 * This class overrides the original implementation to enable loading of ThemeResources (and legacy-styles.css)
 * 
 * 
 * @author klemela
 */
public class ChipsterVaadinServlet extends VaadinServlet {

	@Deprecated
	protected boolean isAllowedVAADINResourceUrl(HttpServletRequest request,
			URL resourceUrl) {
		if ("jar".equals(resourceUrl.getProtocol())) {
			// This branch is used for accessing resources directly from the
			// Vaadin JAR in development environments and in similar cases.

			// Inside a JAR, a ".." would mean a real directory named ".." so
			// using it in paths should just result in the file not being found.
			// However, performing a check in case some servers or class loaders
			// try to normalize the path by collapsing ".." before the class
			// loader sees it.
			
			if (resourceUrl.getPath().contains("!/VAADIN/") || resourceUrl.getPath().contains("!/WebContent/VAADIN/") || resourceUrl.getPath().contains("!/tool-editor/WebContent/VAADIN/")) {
				getLogger().fine(
						"Accepted access to a JAR entry using a class loader: "
								+ resourceUrl);
				return true;
			}
			
			getLogger().info(
					"BLOCKED attempt to access a JAR entry not starting with /VAADIN/: " + resourceUrl + ", getPath(): " + resourceUrl.getPath());
			return false;
			
		} else {
			// Some servers such as GlassFish extract files from JARs (file:)
			// and e.g. JBoss 5+ use protocols vsf: and vfsfile: .

			// Check that the URL is in a VAADIN directory and does not contain
			// "/../"
			if (!resourceUrl.getPath().contains("/VAADIN/")
					|| resourceUrl.getPath().contains("/../")) {
				getLogger().info(
						"Blocked attempt to access file: " + resourceUrl);
				return false;
			}
			getLogger().fine(
					"Accepted access to a file using a class loader: "
							+ resourceUrl);
			return true;
		}
	}

	private static final Logger getLogger() {
		return Logger.getLogger(VaadinServlet.class.getName());
	}
}
