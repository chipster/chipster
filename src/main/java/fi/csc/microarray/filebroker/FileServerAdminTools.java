package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.Strings;

public class FileServerAdminTools {
	public DerbyMetadataServer metadataServer;
	public File cacheRoot;
	public File storageRoot;
	private boolean calculateChecksums;

	public FileServerAdminTools(DerbyMetadataServer metadataServer, File cacheRoot, File storageRoot, boolean calculateChecksums) {
		this.metadataServer = metadataServer;
		this.cacheRoot = cacheRoot;
		this.storageRoot = storageRoot;
		this.calculateChecksums = calculateChecksums;
	}
	
	public String getDataBaseStatusReport() throws SQLException {
		List<String>[] stats = metadataServer.getStatistics();
		
		List<String> keys = stats[0];
		List<String> values = stats[1];
		
		String report = "DATABASE STATUS REPORT\n\n";
		
		for (int i = 0; i < keys.size() && i < values.size(); i++) {
			report += keys.get(i) + "\t" + values.get(i) + "\n";
		}				
		return report;
	}
	
    public String getStorageStatusReport(boolean calculateChecksums) throws UnsupportedEncodingException, FileNotFoundException, FileBrokerException, IOException, SQLException {
    	
    	HashSet<String> dataIds = getAllDataIds();
    	
    	// check which dataIds are ok
    	HashMap<String, String> dataIdErrors = getDataIdErrors(dataIds);
    	HashMap<String, String> dataIdWarnings = getDataIdWarnings(dataIds);
    	    	    	    	    	
    	long errorTotalSize = getTotalSize(dataIdErrors.keySet());
    	long warningTotalSize = getTotalSize(dataIdWarnings.keySet());
    	
    	HashSet<String> dataIdOk = new HashSet<>(dataIds);
    	dataIdOk.removeAll(dataIdErrors.keySet());
    	dataIdOk.removeAll(dataIdWarnings.keySet());
    	
    	long okTotalSize = getTotalSize(dataIdOk);    	    	    
    	
    	HashSet<String> otherFiles = getOtherFiles();
    	long otherFilesTotalSize = getTotalSize(otherFiles);
    	
    	// report
    	
    	String report = "STORAGE STATUS REPORT\n";
    	
    	report += "\n** Summary\n";
    	
    	report += "ok datasets:      \t" + dataIdOk.size() + "\ttotal size: \t" + okTotalSize / 1024 / 1024 + " MB\n";
    	report += "warning datasets: \t" + dataIdWarnings.size() + "\ttotal size: \t" + warningTotalSize / 1024 / 1024 + " MB\n";
    	report += "error datasets:   \t" + dataIdErrors.size() + "\ttotal size: \t" + errorTotalSize / 1024 / 1024 + " MB\n";
    	report += "unknown files:    \t" + otherFiles.size() + "\ttotal size: \t" + otherFilesTotalSize / 1024 / 1024 + " MB\n";
    	
    	if (!dataIdWarnings.isEmpty()) {
    		report += "\n** Warnings\n";
    		for (Entry<String, String> entry : dataIdWarnings.entrySet()) {
    			report += entry.getValue() + "\t" + entry.getKey() + "\n";
    		}
    	}
    	
    	if (!dataIdErrors.isEmpty()) {
    		report += "\n** Errors\n";
    		for (Entry<String, String> entry : dataIdErrors.entrySet()) {
    			report += entry.getValue() + "\t" + entry.getKey() + "\n";
    		}
    	}
    	
    	if (!otherFiles.isEmpty()) {
    		report += "\n** Unknown files\n";
    		for (String filename : otherFiles) {
    			report += filename + "\n";
    		}    	    	   
    	}
    	
    	return report;
	}

    private long getTotalSize(Collection<String> filenames) {
    	
    	long totalSize = 0;    	
    	for (String filename : filenames) {    		
    		File file = getFile(filename);    		
    		totalSize += file.length();
    	}
    	return totalSize;    	
	}

	private HashMap<String, String> getDataIdWarnings(HashSet<String> dataIds) throws UnsupportedEncodingException, FileNotFoundException, IOException {
    	
    	HashMap<String, String> dataIdWarnings = new HashMap<>();
    	  	
    	for (String dataId : dataIds) {    		
    		List<String> warnings = getWarnings(dataId);    		
    		if (!warnings.isEmpty()) {
    			dataIdWarnings.put(dataId, Strings.delimit(warnings, ", "));
    		} 
    	}    	
    	return dataIdWarnings;
	}
    
    private HashMap<String, String> getDataIdErrors(HashSet<String> dataIds) throws UnsupportedEncodingException, FileNotFoundException, IOException, SQLException {
    	
    	HashMap<String, String> dataIdErrors = new HashMap<>();
    	  	
    	for (String dataId : dataIds) {    		
    		List<String> errors = getErrors(dataId);    		
    		if (!errors.isEmpty()) {
    			dataIdErrors.put(dataId, Strings.delimit(errors, ", "));
    		} 
    	}    	
    	return dataIdErrors;
	}

	private List<String> getWarnings(String dataId) throws UnsupportedEncodingException, FileNotFoundException, IOException {

		List<String> warnings = new ArrayList<>();
		
		File file = getFile(dataId);		
		if (file.exists()) {
			String md5 = Md5FileUtils.readMd5WithoutParseException(file);
			if (md5 == null) {
				warnings.add("no md5 file");
			}
		}
		
		return warnings;
    }

	private List<String> getErrors(String dataId) throws UnsupportedEncodingException, FileNotFoundException, IOException, SQLException {
    	    	
    	List<String> errors = new ArrayList<>();

    	File file = getFile(dataId);

    	if (!file.exists()) {
    		errors.add("not on disk");    			
    	}

    	if (calculateChecksums) {
    		String md5 = Md5FileUtils.readMd5WithoutParseException(file);
    		if (md5 != null && file.exists()){
    			if (!md5.equals(Md5FileUtils.calculateMd5(file))) {
    				errors.add("checksum error");
    			}
    		}
    	}

    	DbFile dbFile = metadataServer.fetchFile(dataId);

    	if (dbFile == null) {
    		errors.add("not in database");
    	} else {
    		if (file.exists() && file.length() != dbFile.getSize()) {
    			errors.add("different size");
    		}
    	}

    	return errors;    		    	
	}

	private File getFile(String dataId) {
		return new File(storageRoot, dataId);
	}

	private HashSet<String> getAllDataIds() throws SQLException {
    	// collect dataIds from disk
    	
    	HashSet<String> dataIds = getDataFiles();
    	dataIds.addAll(getMd5Files());
    	    	    	
    	// collect dataIds from db
    	
    	for (DbFile file : metadataServer.listAllFiles()) {
    		dataIds.add(file.getUuid());
    	}
    	
    	return dataIds;
	}

	private HashSet<String> getMd5Files() {		
    	HashSet<String> dataIds = new HashSet<>();    	
    	for (File file : storageRoot.listFiles()) {
    		if (isMd5File(file)) {
    			dataIds.add(file.getName().replace(".md5", ""));
    		}
    	}    
    	return dataIds;
	}

	private HashSet<String> getDataFiles() {
    	HashSet<String> dataIds = new HashSet<>();    	
    	for (File file : storageRoot.listFiles()) {
    		if (isDataFile(file)) {
    			dataIds.add(file.getName());
    		}
    	}    	
    	return dataIds;
	}
	
	private HashSet<String> getOtherFiles() {
    	HashSet<String> filenames = new HashSet<>();    	
    	for (File file : storageRoot.listFiles()) {
    		if (!isDataFile(file) && !isMd5File(file)) {
    			filenames.add(file.getName());
    		}
    	}    	
    	return filenames;
	}
	
	private boolean isDataFile(File file) {		
		return CryptoKey.validateKeySyntax(file.getName());
	}
	
	private boolean isMd5File(File file) {
		return file.getName().endsWith(".md5");
	}
}