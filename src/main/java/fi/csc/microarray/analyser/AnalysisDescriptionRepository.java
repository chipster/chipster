package fi.csc.microarray.analyser;

import java.util.LinkedHashMap;
import java.util.LinkedList;

import org.apache.log4j.Logger;


/**
 * 
 * Any access to descriptions or visibleDescriptions should be use synchronized(this).
 * 
 * @author hupponen
 *
 */
public class AnalysisDescriptionRepository {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(AnalysisDescriptionRepository.class);
	
	private LinkedHashMap<String, AnalysisDescription> descriptions = new LinkedHashMap<String, AnalysisDescription>(); 
	private LinkedHashMap<String, AnalysisDescription> visibleDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedList<AnalysisHandler> handlers = new LinkedList<AnalysisHandler>(); 

	
	
	public void addAnalysisHandler(AnalysisHandler handler) {
		handlers.add(handler);
	}
	
	
	public void loadOperation(String sourceResourceName, boolean hidden) throws AnalysisException {
		AnalysisDescription description = loadDescription(sourceResourceName);
		addDescription(description, hidden);
		
	}
	
	
	private AnalysisDescription loadDescription(String sourceResourceName) throws AnalysisException {
		for (AnalysisHandler handler : handlers) {
			logger.debug("using handler " + handler.getClass().getSimpleName() + ", checking " + sourceResourceName);
			if (handler.supports(sourceResourceName)) {
				return handler.handle(sourceResourceName);
			}
		}
		throw new IllegalArgumentException("none of the loaded handlers support " + sourceResourceName);
	}
	
	
	
	private void addDescription(AnalysisDescription description, boolean hidden) {
		logger.debug("added operation " + description.getFullName());
		
		synchronized(this) {
			descriptions.put(description.getFullName(), description);
			if (!hidden) {
				visibleDescriptions.put(description.getFullName(), description);
			}
		}
	}
	
	public AnalysisDescription getDescription(String fullName) throws AnalysisException {
		AnalysisDescription desc; 

		// get the description
		synchronized(this) {
			desc = descriptions.get(fullName);
		}
		
		// check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			AnalysisDescription newDescription = loadDescription(desc.getSourceResourceName());
			synchronized(this) {
				descriptions.remove(fullName);
				descriptions.put(newDescription.getFullName(), newDescription);
				assert(newDescription.getFullName().equals(fullName));
				if (visibleDescriptions.containsKey(fullName)) {
					visibleDescriptions.remove(fullName);
					visibleDescriptions.put(newDescription.getFullName(), newDescription);
				}
			}
			desc = newDescription;
		}
		
		return desc; 
	}
	
	/**
	 * Returns one huge VVSADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 */
	public StringBuffer serialiseAsStringBuffer() {
		StringBuffer buf = new StringBuffer();
		for (AnalysisDescription description : visibleDescriptions.values()) {
			buf.append(description.getVVSADL());
		}
		return buf;
	}

}
