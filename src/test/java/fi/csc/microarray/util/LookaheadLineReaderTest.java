package fi.csc.microarray.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import org.junit.Assert;
import org.junit.Test;

public class LookaheadLineReaderTest {
	
	@Test
	public void test() throws IOException {
		String buffer = "hello hello\nhello hello\nhello hello\nhowdy\nrockrock";
		LookaheadLineReader lineReader = new LookaheadLineReader(new BufferedReader(new StringReader(buffer)));
		lineReader.peekLine(3);
		lineReader.readLine();
		lineReader.readLine();
		lineReader.read("hello hell".length());
		Assert.assertEquals("o", lineReader.readLine());
		Assert.assertEquals("howdy", lineReader.peekLine());
	}
}
