package fi.csc.microarray.comp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ToolTestHtmlFile {

	private static final String CSS = "<style type=\"text/css\">" + 
			"th {text-align: left; border-bottom-width: 1; border-bottom-style: solid}" +
			"td {padding-right: 1em}" +
			"h3 {margin-top: 2em}" +
			"a {text-decoration: none}" +
			"</style>";

	
	public static void write(ToolTestSummary summary, File webDir) throws IOException {

		File htmlFile = new File(webDir, "index.html");
		File sessionsDir = new File(webDir, SessionReplayTest.SESSIONS_DIR_NAME);
		File screenOutputsDir = new File(webDir, SessionReplayTest.SCREEN_OUTPUTS_DIR_NAME);

		// create files
		FileWriter writer = new FileWriter(htmlFile);
		try {
			writer.write("<html>");
			writer.write("<head>" + CSS + "</head>");

			writer.write("<body>");

			// Main title
			String titleStatus = "<span style=\"color: green\">everything ok!</span>";
			if (summary.getFailedJobs().size() > 0 || summary.getSessionsWithErrors().size() > 0 || summary.getSessionsWithMissingTools().size() > 0) {
				titleStatus = "<span style=\"color: red\">" + 
						summary.getFailCounts().keySet().size() + " tool(s) failed in " + summary.getFailedJobs().size() + " test(s), " +
						summary.getSessionsWithErrors().size() + " session(s) with errors, " +
						summary.getSessionsWithMissingTools().size() + " sessions(s) with missing tools" +
						"</span>";
			}
			writer.write("<h2>Tool tests &ndash; " + titleStatus + "</h2>");

			writer.write("<h3>Summary</h3>");

			// Summary
			
			// times
			String timesHtml = "";
			if (summary.getDuration() > 0) {
				long duration = summary.getDuration()/1000;
				timesHtml += "<tr><td>Total time</td><td>" + String.format("%02dm %02ds", (duration/60), (duration%60)) + "</td></tr>";
			}
			
			if (summary.getStartTime() != null) {
				timesHtml += "<tr><td>Start time</td><td>" + summary.getStartTime() + "</td></tr>";
			}
			writer.write("<table>" +
					"<tr><td>Results summary</td><td>" + 
					summary.getSuccessfulJobs().size() + " <span" + (summary.getSuccessfulJobs().isEmpty() ? "" : " style=\"color: green\"") + ">ok</span>, " + 
					summary.getFailedJobs().size() + " <span" + (summary.getFailedJobs().isEmpty() ? "" : " style=\"color: red\"") + ">failed</span>, "+
					summary.getToolTestResults().size() + " total</td</tr>" +

					"<tr><td>Tool coverage</td><td>" + summary.getUniqueTools().size() + "/" + summary.getAllTools().size() + "</td></tr>" +
					"<tr><td>Sessions</td>" + 
					"<td>" + summary.getUniqueSessions().size() + " total, " + 
					summary.getSessionsWithErrors().size() + " <span" + (summary.getSessionsWithErrors().isEmpty() ? "" : " style=\"color: red\"") + ">with errors</span>, " +
					summary.getSessionsWithMissingTools().size() + " <span" + (summary.getSessionsWithMissingTools().isEmpty() ? "" : " style=\"color: red\"") + ">with missing tools</span> "+
					"</td></tr>" +
					timesHtml +
					"</table>");

			// Sessions with errors
			if (summary.getSessionsWithErrors().size() > 0) {
				writer.write("<h3>Sessions with errors</h3>");
				writer.write("<table><tr>" + 
						"<th>Session</th>" + 
						"<th>Exception</th>" +
						"</tr>");
				for (File f : summary.getSessionsWithErrors().keySet()) {
					writer.write("<tr>" +
							"<td>" + f.getName() + "</td>" +
							"<td>" + "<a href=\"" + SessionReplayTest.STACKTRACES_DIR_NAME + "/" + summary.getSessionsWithErrors().get(f) +"\">" + "exception" + "</a>" + "</td>" +
							"</tr>");
				}
				writer.write("</table>");
			}

			// Missing tools
			if (summary.getSessionsWithMissingTools().size() > 0) {
				writer.write("<h3>Sessions with missing tools</h3>");
				writer.write("<table><tr>" + 
						"<th>Session</th>" + 
						"<th>Tool</th>" +
						"</tr>");
				for (File f : summary.getSessionsWithMissingTools().keySet()) {
					writer.write("<tr>" +
							"<td>" + f.getName() + "</td>" +
							"<td>" + summary.getSessionsWithMissingTools().get(f) + "</td>" +
							"</tr>");
				}
				writer.write("</table>");
			}



			// Tool test results
			writer.write("<h3>Tool test results</h3>");
			writer.write("<table><tr>" + 
					"<th>Tool</th>" + 
					"<th>Tool id</th>" + 
					"<th>Result</th>" +
					"<th>Session</th>" +
					"<th>Task state</th>" + 
					"<th>Test error message</th>" + 
					"<th>Task error message</th>" +
					"<th>Task screen output</th>" + 
					"<th>Duration</th>" + 
					"<th>Outputs with mismatching sizes</th>" + 
					"<th>Outputs with mismatching contents</th>" + 
					"</tr>");

			// Failed tests
			for (ToolTestResult toolTestResult : summary.getFailedJobs()) {
				writer.write("<tr>" +
						"<td>" + toolTestResult.getToolFullName() + "</nobr></td>" +
						"<td><nobr>" + toolTestResult.getToolId() + "</nobr></td>" +
						"<td style=\"color: red\">" + toolTestResult.getTestResult() + "</td>" +
						"<td>" + createSessionLink(toolTestResult.getSession(), sessionsDir) + "</td>" + 
						"<td>" + toolTestResult.getTaskState() + "</td>" +
						"<td>" + nullToEmpty(toolTestResult.getTestErrorMessage()) + "</td>" +
						"<td>" + nullToEmpty(toolTestResult.getTaskErrorMessage()) + "</td>" +
						"<td>" + createScreenOutputLink(toolTestResult.getTaskId(), screenOutputsDir) + "</td>" +
						"<td><nobr>" + getDurationString(toolTestResult.getTaskDuration()) + "</nobr></td>" +
						"<td>" + getOutputsWithMisMatchingSizes(toolTestResult) + "</td>" +
						"<td>" + getOutputsWithMisMatchingContents(toolTestResult) + "</td>" +
						"</tr>");
			}

			// Successful tests
			for (ToolTestResult toolTestResult : summary.getSuccessfulJobs()) {
				writer.write("<tr>" +
						"<td><nobr>" + toolTestResult.getToolFullName() + "</nobr></td>" +
						"<td><nobr>" + toolTestResult.getToolId() + "</nobr></td>" +
						"<td>" + toolTestResult.getTestResult() + "</td>" +
						"<td>" + createSessionLink(toolTestResult.getSession(), sessionsDir) + "</td>" + 
						"<td>" + toolTestResult.getTaskState() + "</td>" +
						"<td>" + nullToEmpty(toolTestResult.getTestErrorMessage()) + "</td>" +
						"<td>" + nullToEmpty(toolTestResult.getTaskErrorMessage()) + "</td>" +
						"<td>" + createScreenOutputLink(toolTestResult.getTaskId(), screenOutputsDir) + "</td>" +
						"<td><nobr>" + getDurationString(toolTestResult.getTaskDuration()) + "</nobr></td>" +
						"<td>" + getOutputsWithMisMatchingSizes(toolTestResult) + "</td>" +
						"<td>" + getOutputsWithMisMatchingContents(toolTestResult) + "</td>" +
						"</tr>");
			}
			writer.write("</table>");


			// Coverage
			writer.write("<h3>Coverage</h3>");
			writer.write("<table><tr>" + 
					"<th>Tool </th>" + 
					"<th>Tool id</th>" + 
					"<th>Count</th>" +
					"<th>Sessions</th>" + 
					"</tr>");


			// Tools with tests
			for (String toolID : summary.getToolsSortedbyTestCount()) {
				String sessionsString = "";
				for (File session : summary.getToolToSessionsMap().get(toolID)) {
					sessionsString += session.getName() + " ";
				}
				sessionsString.trim();

				writer.write("<tr>" +
						//					"<td>" + this.getOperationDefinition(toolID, toolModules).getFullName()  + "</td>" +
						"<td>" + summary.getAllTools().get(toolID) + "</td>" +
						"<td>" + toolID + "</td>" +
						"<td>" + summary.getUniqueTools().get(toolID) + "</td>" +
						"<td>" + sessionsString + "</td>" +
						"</tr>");
			}

			// Tools without tests
			for (String toolId : summary.getAllTools().keySet()) {
				if (!summary.getUniqueTools().containsKey(toolId)) {
					writer.write("<tr>" +
							"<td>" + summary.getAllTools().get(toolId) + "</td>" +
							"<td>" + toolId + "</td>" +
							"<td>" + 0 + "</td>" +
							"<td>" + "" + "</td>" +
							"</tr>");
				}
			}


			writer.write("</table>");


			writer.write("</body></html>");
			writer.flush();
		} finally {
			writer.close();
		}
	}


	
	private static String nullToEmpty(String s) {
		if (s == null) {
			return "";
		} else {
			return s;
		}
	}
	
	private static String createSessionLink(File sessionFile, File sessionsDir) {
		// TODO
		//		return "<a href=\"" + sessionsDir + "/" + sessionFile.getName() + "\">" + sessionFile.getName() + "</a>";

		return sessionFile.getName();
	}
	
	
	private static String createScreenOutputLink(String taskId, File screenOutputsDir) {
		return "<a href=\"" + screenOutputsDir.getName() + "/" + taskId  + SessionReplayTest.SCREEN_OUTPUT_POSTFIX +"\">" + "output" + "</a>";
	}

	private static String getDurationString(long duration) {
		long d = duration / 1000;
		return String.format("%02dm %02ds", (d/60), (d%60));
	}


	private static String getOutputsWithMisMatchingSizes(ToolTestResult toolTestResult) {
		String s = "";
		for (String output : toolTestResult.getOutputsWithMisMatchingSizes()) {
			s += output + ", ";
		}
		if (s.endsWith(", ")) {
			s = s.substring(0, s.length() - ", ".length());
		}
		return s;
	}

	private static String getOutputsWithMisMatchingContents(ToolTestResult toolTestResult) {
		String s = "";
		for (String output : toolTestResult.getOutputsWithMisMatchingContents()) {
			s += output + ", ";
		}
		if (s.endsWith(", ")) {
			s = s.substring(0, s.length() - ", ".length());
		}
		return s;
	}

	

}

	
