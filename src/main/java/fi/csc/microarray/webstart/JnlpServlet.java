package fi.csc.microarray.webstart;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.XmlUtil;

/**	
 * Servlet for generating modified jnlp files according to url query parameters.
 * 
 * Don't use this with sensitive information, because javaws caches all jnlp files. Current 
 * Java version (1.7) seems to delete application shortcuts when there is a question 
 * mark in the href attribute in the jnlp file.
 * 
 * @author klemela
 */
public class JnlpServlet extends HttpServlet {
	
	@Override
	protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws javax.servlet.ServletException, java.io.IOException {
		
		resp.setContentType("application/x-java-jnlp-file");
		
		String memoryString = req.getParameter("memory");
		
		Integer memory = null;
		if (memoryString != null) {
			try {
				memory = Integer.parseInt(memoryString);
			} catch (NumberFormatException e) {
				resp.sendError(HttpServletResponse.SC_BAD_REQUEST, "parameter 'memory' must be an integer");
				return;
			}
		}
		
		resp.setStatus(HttpServletResponse.SC_OK);
		
		PrintWriter writer = resp.getWriter();
		
		Document jnlp;
		try {
			jnlp = getJnlp(memory);
			XmlUtil.printXml(jnlp, writer);		
		} catch (Exception e) {
			resp.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR, "jnlp modification failed");
			return;
		}
		
		
	}

	private Document getJnlp(Integer memory) throws SAXException, IOException, ParserConfigurationException {
		
		File wsConfigFile = new File(DirectoryLayout.WEB_ROOT + File.separator + "chipster.jnlp");
		if (!wsConfigFile.exists()) {
			throw new IOException("chipster.jnlp not found");
		}
		
		Document doc = XmlUtil.parseFile(wsConfigFile);		
		Element jnlp = (Element)doc.getDocumentElement();
		
		/* When javaws installs a application, it doesn't use the jnlp file that was just downloaded, but
		 * downloads the jnlp again from the url denoted by the file's href attribute. We must generate a query string 
		 * that calls this servlet with the exactly same parameters again. This way the parameters are relayed 
		 * reliably regardless whether the application was already installed or not.
		 */
		String href="servlet.jnlp";
		if (memory != null) {
			href += "?";
		}
		if (memory != null) {
			href += "memory=" + memory;
		}

		jnlp.setAttribute("href", href);

		if (memory != null) {
			Element resources = (Element) jnlp.getElementsByTagName("resources").item(0);
			Element j2se = (Element) resources.getElementsByTagName("j2se").item(0);
			j2se.setAttribute("java-vm-args", "-Xmx" + memory + "m");
		}
		
//		// use this to give parameters for main method
//		Element applicationDesc = (Element)jnlp.getElementsByTagName("application-desc").item(0);
//		
//		Element keyArgument = doc.createElement("argument");
//		Element valueArgument = doc.createElement("argument");
//		keyArgument.setTextContent("-parameter-name");
//		valueArgument.setTextContent(parameterValue);		
//		applicationDesc.appendChild(keyArgument);
//		applicationDesc.appendChild(valueArgument);
		
		return doc;
	}
}
