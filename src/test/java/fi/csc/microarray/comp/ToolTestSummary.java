package fi.csc.microarray.comp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.ParseException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import fi.csc.microarray.comp.SessionReplayTest.TestResult;

public class ToolTestSummary {

	private static final String SESSIONS_WITH_ERRORS_FILE = "sessions-with-errors.txt";
	private static final String SESSIONS_WITH_MISSING_TOOLS_FILE = "sessions-with-missing-tools.txt";
	private static final String ALL_TOOLS_FILE = "all-tools.txt";
	public static final String STARTTIME_FILE_NAME = "starttime.txt";
	public static final String ENDTIME_FILE_NAME = "endtime.txt";
	
	private List<ToolTestResult> toolTestResults = new LinkedList<ToolTestResult>();
	
	// Sessions which cause something to be thrown, normal tool failures etc not counted
	private LinkedHashMap<File, String> sessionsWithErrors = new LinkedHashMap<File, String>();
	
	private LinkedHashMap<File, String> sessionsWithMissingTools = new LinkedHashMap<File, String>();
	final HashMap<String, Integer> uniqueTools = new HashMap<String,Integer>(); 
	private HashMap<String, Integer> failCounts = new HashMap<String,Integer>(); 
	private List<ToolTestResult> failedJobs = new LinkedList<ToolTestResult>();
	private List<ToolTestResult> successfulJobs = new LinkedList<ToolTestResult>();
	private Set<File> uniqueSessions = new HashSet<File>();
	private String[] toolsSortedbyTestCount;
	private HashMap<String, List<File>> toolToSessionsMap = new HashMap<String, List<File>>();
	private LinkedHashMap<String, String> allTools = new LinkedHashMap<String, String>();
 
	private Date startTime;
	private Date endTime;
	
	
	public void calculateStats() {
		countSuccessfulAndFailedJobs();
		
		// count unique tools and sessions

		for (ToolTestResult toolTestResult : toolTestResults) {
			if (!uniqueTools.containsKey(toolTestResult.getToolId())) {
				uniqueTools.put(toolTestResult.getToolId(), 1);
			} else {
				uniqueTools.put(toolTestResult.getToolId(), uniqueTools.get(toolTestResult.getToolId()) + 1);
			}

			if (!toolToSessionsMap.containsKey(toolTestResult.getToolId())) {
				List<File> sessionsList = new LinkedList<File>();
				sessionsList.add(toolTestResult.getSession());
				toolToSessionsMap.put(toolTestResult.getToolId(), sessionsList);
			} else {
				List<File> sessionsList = toolToSessionsMap.get(toolTestResult.getToolId());
				if (!sessionsList.contains(toolTestResult.getSession())) {
					sessionsList.add(toolTestResult.getSession());
				}
			}

			if (!uniqueSessions.contains(toolTestResult.getSession())) {
				uniqueSessions.add(toolTestResult.getSession());
			}
		}

		// sort tools by test count
		toolsSortedbyTestCount = uniqueTools.keySet().toArray(new String[]{});
		Arrays.sort(toolsSortedbyTestCount, new Comparator<String>() {
			@Override
			public int compare(String arg0, String arg1) {
				return uniqueTools.get(arg1) - uniqueTools.get(arg0);
			}
		});

	}
	
	private void countSuccessfulAndFailedJobs() {

		// failed and successful tasks
		for (ToolTestResult toolTestResult : toolTestResults) {
			if (TestResult.FAIL.equals(toolTestResult.getTestResult())) {
				failedJobs.add(toolTestResult);
				String toolID = toolTestResult.getToolId();
				if (failCounts.containsKey(toolID)) {
					failCounts.put(toolID, failCounts.get(toolID) + 1);
				} else {
					failCounts.put(toolID, 1);
				}
			} else {
				successfulJobs.add(toolTestResult);
			}
		}
	}

	public List<ToolTestResult> getToolTestResults() {
		return this.toolTestResults;
	}

	public List<ToolTestResult> getFailedJobs() {
		return this.failedJobs;
	}
	
	public List<ToolTestResult> getSuccessfulJobs() {
		return this.successfulJobs;
	}


	private void writeSessionsWithErrors(File webDir) throws IOException {
		File results = new File(webDir, SESSIONS_WITH_ERRORS_FILE);
		Writer writer = new BufferedWriter(new FileWriter(results));
		try {
			for (File f : sessionsWithErrors.keySet()) {
				writer.write(f.getName() + ToolTestTextFile.TEXT_FILE_SEPARATOR + sessionsWithErrors.get(f) + "\n");
			}
		} finally {
			writer.close();
		}
	}

	private void writeSessionsWithMissingTools(File webDir) throws IOException {
		File results = new File(webDir, SESSIONS_WITH_MISSING_TOOLS_FILE);
		Writer writer = new BufferedWriter(new FileWriter(results));
		try {
			for (File f : sessionsWithMissingTools.keySet()) {
				writer.write(f.getName() + ToolTestTextFile.TEXT_FILE_SEPARATOR + sessionsWithMissingTools.get(f) + "\n");
			}
		} finally {
			writer.close();
		}
	}

	private void writeAllTools(File webDir) throws IOException {
		File f = new File(webDir, ALL_TOOLS_FILE);
		Writer writer = new BufferedWriter(new FileWriter(f));
		try {
			for (String toolId : allTools.keySet()) {
				writer.write(toolId + ToolTestTextFile.TEXT_FILE_SEPARATOR + allTools.get(toolId) + "\n");
			}
		} finally {
			writer.close();
		}
	}

	private void writeTimes(File webDir) throws IOException {
		BufferedWriter writer;
		
		if (startTime != null) {
			writer = new BufferedWriter(new FileWriter(new File(webDir, STARTTIME_FILE_NAME)));
			try {
				writer.write(Long.toString(startTime.getTime()));
			} finally {
				writer.close();
			}
		}
		
		if (endTime != null) {
			writer = new BufferedWriter(new FileWriter(new File(webDir, ENDTIME_FILE_NAME)));
			try {
				writer.write(Long.toString(endTime.getTime()));
			} finally {
				writer.close();
			}
		}
	}
		
	
	private void readAllTools(File webDir) throws IOException {
		allTools = new LinkedHashMap<String, String>();
		File f = new File(webDir, ALL_TOOLS_FILE);
		if (!f.exists()) {
			return;
		}

		BufferedReader reader = new BufferedReader(new FileReader(f));
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().isEmpty()) {
					continue;
				}

				String[] parts = line.split(ToolTestTextFile.TEXT_FILE_SEPARATOR, -1);
				allTools.put(parts[0], parts[1]);		
			}
		} finally {
			reader.close();
		}
	}

	
	private void readTimes(File webDir) throws IOException, ParseException {
		File f = new File(webDir, STARTTIME_FILE_NAME);
		if (!f.exists()) {
			return;
		}

		BufferedReader reader = new BufferedReader(new FileReader(f));
		try {
			String line = reader.readLine();
			if (!line.trim().isEmpty()) {
				startTime = new Date(Long.parseLong(line));
			}
		} finally {
			reader.close();
		}

		f = new File(webDir, ENDTIME_FILE_NAME);
		if (!f.exists()) {
			return;
		}

		reader = new BufferedReader(new FileReader(f));
		try {
			String line = reader.readLine();
			if (!line.trim().isEmpty()) {
				endTime = new Date(Long.parseLong(line));
			}
		} finally {
			reader.close();
		}
	}
	
	private void readSessionsWithErrors(File webDir) throws IOException {
		sessionsWithErrors = new LinkedHashMap<File, String>();
		
		File results = new File(webDir, SESSIONS_WITH_ERRORS_FILE);
		if (!results.exists()) {
			return;
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(results));
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().isEmpty()) {
					continue;
				}
				String[] parts = line.split(ToolTestTextFile.TEXT_FILE_SEPARATOR, -1);
				sessionsWithErrors.put(new File(parts[0]), parts[1]);
			}
		} finally {
			reader.close();
		}
	}

	private void readSessionsWithMissingTools(File webDir) throws IOException {
		sessionsWithMissingTools = new LinkedHashMap<File, String>();

		File results = new File(webDir, SESSIONS_WITH_MISSING_TOOLS_FILE);
		if (!results.exists()) {
			return;
		}

		BufferedReader reader = new BufferedReader(new FileReader(results));
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.trim().isEmpty()) {
					continue;
				}
				String[] parts = line.split(ToolTestTextFile.TEXT_FILE_SEPARATOR, -1);
				sessionsWithMissingTools.put(new File(parts[0]), parts[1]);
			}
		} finally {
			reader.close();
		}
	}

	public void writeToFiles(File webDir) throws IOException {
		ToolTestTextFile.write(this, webDir);
		ToolTestHtmlFile.write(this, webDir);
		writeSessionsWithErrors(webDir);
		writeSessionsWithMissingTools(webDir);
		writeAllTools(webDir);
		writeTimes(webDir);
	}

	public void readFromFiles(File webDir) throws IOException, ParseException {
		this.toolTestResults = ToolTestTextFile.parse(webDir);
		readSessionsWithErrors(webDir);
		readSessionsWithMissingTools(webDir);
		readAllTools(webDir);
		readTimes(webDir);
	}

	public HashMap<String, Integer> getUniqueTools() {
		return uniqueTools;
	}

	public HashMap<String, Integer> getFailCounts() {
		return failCounts;
	}

	public LinkedHashMap<File, String> getSessionsWithErrors() {
		return sessionsWithErrors;
	}

	public LinkedHashMap<File, String> getSessionsWithMissingTools() {
		return sessionsWithMissingTools;
	}

	public Set<File> getUniqueSessions() {
		return uniqueSessions;
	}

	public String[] getToolsSortedbyTestCount() {
		return toolsSortedbyTestCount;
	}

	public HashMap<String, List<File>> getToolToSessionsMap() {
		return toolToSessionsMap;
	}

	public LinkedHashMap<String, String> getAllTools() {
		return allTools;
	}

	public Date getStartTime() {
		return startTime;
	}

	public void setStartTime(Date startTime) {
		this.startTime = startTime;
	}

	public Date getEndTime() {
		return endTime;
	}

	public void setEndTime(Date endTime) {
		this.endTime = endTime;
	}

	public long getDuration() {
		if (startTime != null && endTime != null) {
			return endTime.getTime() - startTime.getTime();
		} else {
			return -1;
		}
	}
}

