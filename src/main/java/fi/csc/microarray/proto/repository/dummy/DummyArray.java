package fi.csc.microarray.proto.repository.dummy;

import java.io.ByteArrayInputStream;
import java.io.InputStream;

import fi.csc.microarray.proto.repository.Array;

public class DummyArray implements Array {

	private String name;
	
	public DummyArray(String name) {
		this.name = name;
	}
		
	public Platform getPlatform() {
		return Platform.CDNA;
	}

	public String getName() {
		return name;
	}

	public InputStream getContents() {
		return new ByteArrayInputStream("contents".getBytes());
	}

}
