package fi.csc.microarray.proto.repository;

import java.io.InputStream;

public interface Array {
	
	public static enum Platform {
		CDNA,
		AFFYMETRIX,
		OTHER;
	}
	
	public String getName();
	
	public Platform getPlatform();
	
	public InputStream getContents();

}
