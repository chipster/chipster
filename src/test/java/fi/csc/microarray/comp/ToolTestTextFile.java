package fi.csc.microarray.comp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.comp.SessionReplayTest.TestResult;

public class ToolTestTextFile {

	public static final String TEXT_FILE_SEPARATOR = ";;;;";
	private static final String TEXT_RESULTS_FILE = "results.txt";
	
	public static void write(ToolTestSummary summary, File webDir) throws IOException {
		
		File resultsFile = new File(webDir, TEXT_RESULTS_FILE);
		FileWriter writer = new FileWriter(resultsFile);
		try {
			for (ToolTestResult result : summary.getToolTestResults()) {
				writer.write(
						result.getTestResult() + TEXT_FILE_SEPARATOR +
						result.getToolName() + TEXT_FILE_SEPARATOR +
						result.getToolId() + TEXT_FILE_SEPARATOR +
						result.getToolFullName() + TEXT_FILE_SEPARATOR +
						result.getTaskId() + TEXT_FILE_SEPARATOR +
						result.getTaskState() + TEXT_FILE_SEPARATOR +
						nullToEmpty(result.getTaskErrorMessage()).replaceAll("\n", " ") + TEXT_FILE_SEPARATOR +
						result.getTaskDuration() + TEXT_FILE_SEPARATOR +
						result.getSession() + TEXT_FILE_SEPARATOR +
						nullToEmpty(result.getTestErrorMessage()) + TEXT_FILE_SEPARATOR
						);
				writer.write("\n");
			}		
		} finally {
			org.apache.commons.io.IOUtils.closeQuietly(writer);
		}
	}

	public static List<ToolTestResult> parse(File inputDir) throws IOException {
		List<ToolTestResult> toolTestResults = new LinkedList<ToolTestResult>();
		File resultsFile = new File(inputDir, TEXT_RESULTS_FILE);
		if (!resultsFile.exists()) {
			throw new FileNotFoundException(resultsFile.getAbsolutePath());
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(resultsFile));
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().isEmpty()) {
					continue;
				}
				
				String[] parts = line.split(TEXT_FILE_SEPARATOR, -1);
				TestResult testResult = TestResult.valueOf(parts[0]); 
				String toolName = parts[1]; 
				String toolId = parts[2]; 
				String toolFullName = parts[3]; 
				String taskId = parts[4]; 
				String taskState = parts[5]; 
				String taskErrorMessage = parts[6]; 
				long taskDuration = Long.parseLong(parts[7]); 
				File sessionFile = new File(parts[8]); 
				String testErrorMessage = parts[9]; 
				
				ToolTestResult singleResult = new ToolTestResult(testResult, taskId, taskState, taskErrorMessage, taskDuration, 
															toolId, toolName, toolFullName, 
															sessionFile, testErrorMessage);
				toolTestResults.add(singleResult);
				
			}
		} finally {
			org.apache.commons.io.IOUtils.closeQuietly(reader);
		}
		return toolTestResults;
	}
	
	private static String nullToEmpty(String s) {
		if (s == null) {
			return "";
		} else {
			return s;
		}
	}

	
}
