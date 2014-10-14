package fi.csc.microarray.client.dataimport;

import java.io.File;
import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.util.SystemMonitorUtil;

public class DataParsingTest {
	
	private static final int MB = 1024*1024;
	
	class TestInformator implements ProgressInformator{

		private long memoryLimit;
		
		public TestInformator(long memoryLimit) {
			this.memoryLimit = memoryLimit;
		}
		
		private int max;
		
		public void destroyInformator() { }

		public void initializeInformator() { }

		public void setMaximumValue(int max) {
			this.max = max;
		}

		public void setMessage(String message) {
			System.out.println("Job informator: " + message);
		}

		public void setMinimunValue(int min) { }

		public void setProcess(RunnableImportProcess process) { }

		public void setValue(int state) {
			// Every 10000
			if(state % 10000 == 0 && state != 0){
				System.out.println("Line: " + state + "   percents done: " + (int)(((float)state/(float)max) * 100) + "%   memory usage: " + SystemMonitorUtil.getMemInfo());
			}
			
			// If memory usage raises too high the test is failed
			Assert.assertTrue("Memory limit exceeded. Limit: " + SystemMonitorUtil.bytesToMegas(memoryLimit) + " Mb, used " + SystemMonitorUtil.bytesToMegas(SystemMonitorUtil.getUsed()), SystemMonitorUtil.getUsed() < memoryLimit);
			
		}

		public void stopProcess() { }

		public void setIndeterminate(boolean newValue) {
			// Do nothing
		}
	}
	
	@Test
	public void testSmallSizedRealData() throws IOException{
		System.out.println("\n");
		System.out.println("Small sized real data (Rows: ~15000, Column: 4, Size: 384kb, Delimiter: Tab, Decimal separator: Dot)");
		System.out.println("Showing all rows and columns");
		System.out.println("===============================");
		long start = System.currentTimeMillis();
		ConversionModel model = new ConversionModel(null);
		model.setInputFile(new File("src/test/resources/affy_example.cel"));
		//originally the limit was 8 MB so probably it has been enough at some point, but now we need more
		model.chopData(false, new TestInformator(30*MB));
		System.out.println("Best suitable delimiter was: " + model.getDelim().getName());
		System.out.println("===============================");
		System.out.println("Parsing done. Total time: " + (System.currentTimeMillis() - start));
	}
}
