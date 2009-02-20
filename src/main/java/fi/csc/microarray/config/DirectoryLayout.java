package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;


/**
 * Specifies Chipster directory layout. This class if the decisive specification 
 * for directory layout. Directory layout is derived using three sources of 
 * information: 1) constant paths in this class, 2) detected runtime platform,
 * and 3) configuration. For a human readable description of directory layout see 
 * <a href="http://chipster.wiki.sourceforge.net/DirectoryLayout">http://chipster.wiki.sourceforge.net/DirectoryLayout</a>. 
 * 
 * @author Aleksi Kallio
 *
 */
public class DirectoryLayout {

	public enum Type {
		CLIENT,
		SERVER;
	}

	private Type type;
	
	public DirectoryLayout(Type type) {
		this.type = type;
	}

	public File getFileroot() throws IOException {
		if (type == Type.SERVER) {
			File fileRepository = new File(MicroarrayConfiguration.getValue("frontend", "fileServerPath"));
			if (!fileRepository.exists()) {
				boolean ok = fileRepository.mkdir();
				if (!ok) {
					throw new IOException("could not create file root at " + fileRepository);
				}
			}
			return fileRepository;
			
		} else {
			throw new UnsupportedOperationException();
		}
	}
	
}
