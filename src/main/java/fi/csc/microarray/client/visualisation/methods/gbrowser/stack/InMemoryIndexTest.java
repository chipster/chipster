package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.File;

import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GtfUtil;

public class InMemoryIndexTest {
	
	public static void main(String[] args) {
		long t = System.currentTimeMillis();
		System.out.println("Start mem: " + getMemoryUsage());
		File file = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66.gtf");
		GtfUtil.loadFile(file);
		
		System.out.println("End mem: " + getMemoryUsage());
		System.out.println("Time: " + (System.currentTimeMillis() - t) / 1000 + " s.");
	}
	
	private static long getMemoryUsage() {
		Runtime rt = Runtime.getRuntime();
		return (rt.totalMemory() - rt.freeMemory()) / 1024 / 1024;
	}
}