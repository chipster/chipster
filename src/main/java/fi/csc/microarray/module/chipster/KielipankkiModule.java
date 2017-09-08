package fi.csc.microarray.module.chipster;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontFormatException;
import java.awt.GraphicsEnvironment;
import java.io.IOException;
import java.util.Collections;
import java.util.Enumeration;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import javax.jms.JMSException;
import javax.swing.Icon;
import javax.swing.JMenu;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;

/**
 * Client module for The Language Bank of Finland called Kielipankki
 * https://www.kielipankki.fi
 * 
 * @author klemela
 *
 */
public class KielipankkiModule implements Module {
	
	private static Logger logger = LogManager.getLogger(KielipankkiModule.class);
	
	public static class TypeTags {
	}
	
	public static class VisualisationMethods {
	}

	public void plugContentTypes(DataManager manager) {
	}
	

	public void plugFeatures(DataManager manager) {
	}

	public void plugModifiers(DataManager manager) {
		// nothing to plug
	}

	@Override
	public String[] getServerModuleNames() {
		return new String[] { "kielipankki" };
	}

	@Override
	public String getModuleLongName(String moduleName) {
		if ("kielipankki".equals(moduleName)) {
			return "Kielipankki";
		} else {
			return moduleName;
		}
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
	public boolean isImportToolSupported() {
		return false;
	}

	@Override
	public boolean isWorkflowCompatible(DataBean data) {
		// TODO replace inclusive check with exclusive: check for raw expression data and other illegal types
//		return ChipsterInputTypes.GENE_EXPRS.isTypeOf(data);
		return true;
	}

	@Override
	public VisualisationMethod[] getVisualisationMethods() {
		return new VisualisationMethod[] {
		};
	}

	@Override
	public List<DbSession> getExampleSessions(boolean isStandalone) throws JMSException, Exception {
		
		 List<DbSession> sessions = Session.getSession().getServiceAccessor().getFileBrokerClient().listPublicRemoteSessions();
		 
		return sessions;
	}

	@Override
	public String[][] getRepositoryWorkflows() {
		return new String[][] { 
		};
	}

	@Override
	public boolean isMetadata(DataBean data) {
		return false;
	}
	
	@Override
	public void preProcessInputMetadata(Operation oper, DataBean metadataInput) throws MicroarrayException, IOException {
	}

	@Override
	public void postProcessOutputMetadata(OperationRecord oper, DataBean metadataOutput) throws MicroarrayException, IOException {
	}

	@Override
	public String getShortDataName(String categoryName) {
		return BasicModule.shortenDataName(categoryName);
	}

	@Override
	public boolean countOperationResults() {
		return true;
	}

	/**
	 * Generates nice context link panel for quickly using genome browser. If not in standalone
	 * mode, null is returned. 
	 */
	@SuppressWarnings("serial")
	@Override
	public JPanel getContextLinkPanel(int selectedDataCount) {

		return null;
	}

	@Override
	public boolean notesVisibleAtStartup() {
		return true;
	}

	@Override
	public String getDisplayName() {
		return "Mylly";
	}

	@Override
	public String getManualHome() {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		return configuration.getString("client", "manual-root");
	}

	@Override
	public List<Boolean> flagLinkableColumns(Table columns, DataBean data) {
		return Collections.nCopies(columns.getColumnCount(), false);
	}

	@Override
	public IntegratedEntity createLinkableEntity(Table columns, DataBean data) {
		return null;
	}

	@Override
	public void addTypeTags(DataBean data) throws MicroarrayException, IOException {
	}

	@Override
	public Icon getIconFor(DataBean data) {
		return data.getContentType().getIcon();
	}


	@Override
	public void updateUIDefaults(UIDefaults defaults) {
		defaults.put("SimpleInternalFrame.activeTitleBackground", new Color(77, 90, 145));
		defaults.put("ScrollBar.thumb", new Color(122, 144, 195));
		
		// replace all l&f settings with the current background color with the Kielipankki background color
		Color defaultBackground = UIManager.getColor("Panel.background");
		Color background = new Color(202, 210, 230);
		
    	for (Enumeration<Object> e = UIManager.getDefaults().keys(); e.hasMoreElements(); ){
    		Object obj = e.nextElement();
    		if (obj instanceof String){
    			String key = (String)obj;
    			Color value = UIManager.getColor(key);
    			if (defaultBackground.equals(value)) {
    				defaults.put(key, background);
    			}
    		}
    	}
    	
    	// change font (breaks text size changes)
    	try {
			GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
			Font font = Font.createFont(Font.TRUETYPE_FONT, KielipankkiModule.class.getResourceAsStream("/fonts/Lato-Regular.ttf"));
			font = font.deriveFont(((SwingClientApplication)Session.getSession().getApplication()).getFontSize() + 2);
			ge.registerFont(font);
			
			for (Entry<Object, Object> entry : UIManager.getDefaults().entrySet()) {
				Object key = entry.getKey();
				Object value = javax.swing.UIManager.get(key);
				if (value != null && value instanceof javax.swing.plaf.FontUIResource) {
					// Lato font doesn't render OSX symbols correctly
					if (!"MenuItem.acceleratorFont".equals(key.toString())) {
						javax.swing.UIManager.put(key, new FontUIResource(font));
					}	
				}
			}
		} catch (IOException|FontFormatException e) {
			// creating an error dialog would fail
			logger.error("font loading error", e);
			System.err.println("font loading error: " + e);
		}
	}
}
