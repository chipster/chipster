package fi.csc.microarray.comp;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.comp.SessionReplayTest.TestResult;

class ToolTestResult {
	private TestResult testResult;
	private String taskId;
	private String taskState;
	private String taskErrorMessage;
	private long taskDuration;
	
	
	
	private String toolId;
	private String toolName;
	private String toolFullName;
	private File sessionFile;
	private String testErrorMessage;
	private List<String> outputsWithMisMatchingSizes = new LinkedList<String>();
	private List<String> outputsWithMisMatchingContents = new LinkedList<String>();


	public ToolTestResult(TestResult testResult, String taskId, String taskState, String taskErrorMessage, long taskDuration, 
			String toolId, String toolName, String toolFullName, File sessionFile, String testErrorMessage) {
		this.testResult = testResult;
		this.taskId = taskId;
		this.taskState = taskState;
		this.taskErrorMessage = taskErrorMessage;
		this.taskDuration = taskDuration;
		
		this.toolId = toolId;
		this.toolName = toolName;
		this.toolFullName = toolFullName;
		
		this.sessionFile = sessionFile;
		this.testErrorMessage = testErrorMessage;
	}

	
	public ToolTestResult(TestResult testResult, File sessionFile, Task task, String testErrorMessage) {
		this.testResult = testResult;
		this.taskId = task.getId();
		this.taskState = task.getState().name();
		this.taskErrorMessage = task.getErrorMessage();
		this.taskDuration = task.getExecutionTime();
		
		this.toolId = task.getOperationID();
		this.toolName = task.getName();
		this.toolFullName = task.getFullName();
		
		this.sessionFile = sessionFile;
		this.testErrorMessage = testErrorMessage;
	}


	public long getTaskDuration() {
		return taskDuration;
	}

	public String getTaskErrorMessage() {
		return taskErrorMessage;
	}

	
	public String getTaskState() {
		return taskState;
	}

	
	public String getTaskId() {
		return taskId;
	}

	public String getToolId() {
		return toolId;
	}

	public String getToolName() {
		return toolName;
	}

	public String getToolFullName() {
		return toolFullName;
	}

	public TestResult getTestResult() {
		return testResult;
	}
	public File getSession() {
		return sessionFile;
	}
	
	public String getTestErrorMessage() {
		return testErrorMessage;
	}

	public List<String> getOutputsWithMisMatchingSizes() {
		return outputsWithMisMatchingSizes;
	}

	public void setOutputsWithMisMatchingSizes(List<String> outputs) {
		this.outputsWithMisMatchingSizes = outputs;
	}

	public List<String> getOutputsWithMisMatchingContents() {
		return outputsWithMisMatchingContents;
	}

	public void setOutputsWithMisMatchingContents(List<String> outputs) {
		this.outputsWithMisMatchingContents = outputs;
	}

	
}