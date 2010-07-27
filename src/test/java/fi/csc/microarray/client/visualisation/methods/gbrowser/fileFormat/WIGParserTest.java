package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class WIGParserTest {
	
	String path = "src/test/resources/wig";
	
	@Test
	public void testFixedStep() throws IOException {
		
		File wigFile = new File(path, "fixedStep.wig");
		WIGParser wig = new WIGParser(wigFile);
		
		Assert.assertEquals(wig.getHeaderLength(new File(path, "fixedStep.wig")), 86, "incorrect header length");
		
		ChunkDataSource data = new ChunkDataSource(wigFile, wig);
		byte[] fileChunk = new byte[(int)wigFile.length()];
		data.read(wig.getHeaderLength(wigFile), fileChunk);
		List<ColumnType> columns = Arrays.asList(new ColumnType[] {
				ColumnType.VALUE});
		
		//wigFile = new File(path, "fixedStep.wig");
		List<RegionContent> abc = data.getFileParser().getAll(new String(fileChunk), columns);
		Assert.assertEquals(abc.get(0).region.end.bp, new Long(200));
		Assert.assertEquals(abc.get(6).values.get(ColumnType.VALUE).toString(), "1.1013745");
		Assert.assertEquals(wig.getBpRegion(getData(wigFile)).end.bp, new Long(3400));
		Assert.assertEquals(abc.get(2).region.start.bp, new Long(401));
		Assert.assertEquals(abc.get(140).region.end.chr.toString(), "21");
		
	}
	
	public String getData(File file) throws IOException {
		
		FileReader fReader = new FileReader(file);
		BufferedReader reader = new BufferedReader(fReader);
		
		String string = "";
		reader.readLine();
		reader.readLine();
		String tmp = reader.readLine();
		
		while (tmp != null) {
			
			string += tmp+"\n";
			tmp = reader.readLine();
			
		}
		
		return string;
	}
	
	public void testVariableStep() throws IOException {
		
		File wigFile = new File(path, "variableStep.wig");
		WIGParser wig = new WIGParser(wigFile);
		
		Assert.assertEquals(wig.getHeaderLength(new File(path, "variableStep.wig")), 175, "incorrect header length");
		
		ChunkDataSource data = new ChunkDataSource(wigFile, wig);
		byte[] fileChunk = new byte[(int)wigFile.length()];
		data.read(wig.getHeaderLength(wigFile), fileChunk);
		List<ColumnType> columns = Arrays.asList(new ColumnType[] {
				ColumnType.VALUE});
		
		//wigFile = new File(path, "variableStep.wig");
		List<RegionContent> abc = data.getFileParser().getAll(new String(fileChunk), columns);
		Assert.assertEquals(abc.get(0).region.end.bp, new Long(830099));
		Assert.assertEquals(abc.get(6).values.get(ColumnType.VALUE).toString(), "1.27");
		Assert.assertEquals(wig.getBpRegion(getData(wigFile)).end.bp, new Long(851424));
		Assert.assertEquals(abc.get(2).region.start.bp, new Long(830125));
		Assert.assertEquals(abc.get(24).region.end.chr.toString(), "2");
		
	}
	
	public static void main(String[] args) throws Exception {
		WIGParserTest test = new WIGParserTest();
		test.testFixedStep();
		test.testVariableStep();
	}

}