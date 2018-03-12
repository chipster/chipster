package fi.csc.microarray.comp;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.TimeUnit;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;
import io.reactivex.subjects.PublishSubject;


public class ProcessMonitor implements Runnable {

	private static final int MAX_SCREEN_OUTPUT_SIZE = 100000; // number of chars
	private static final int MAX_SCREEN_OUTPUT_UPDATE_SIZE = 10000; // number of chars

	private static final String CLIP_CLIP_LINE = 
			"-----clip-----clip-----clip-----clip-----clip-----clip-----clip-----clip-----\n";

	private Process process;
	private Consumer<String> updateScreenOutputCallback;
	private BiConsumer<JobState, String> finishCallback;
	private Pattern successStringPattern;

	private PublishSubject<Boolean> screenOutputSubject = PublishSubject.create();

	private StringBuffer screenOutput = new StringBuffer();
	private boolean outputBufferFull = false;
	private boolean outputUpdateBufferFull = false;

	static final Logger logger = Logger.getLogger(ProcessMonitor.class);

	public ProcessMonitor(
			Process process, 
			Consumer<String> updateScreenOutputCallback,
			BiConsumer<JobState, String> finishCallback,
			Pattern successStringPattern) {
		this.process = process;
		this.updateScreenOutputCallback = updateScreenOutputCallback;
		this.finishCallback = finishCallback;
		this.successStringPattern = successStringPattern;
	}

	public void run() {

		// throttle screen output updates
		screenOutputSubject
		.throttleLast(1, TimeUnit.SECONDS)
		.subscribe(arg -> {
			if (outputUpdateBufferFull) {
				return;
			}

			String s = screenOutput.toString();
			if (s.length() > MAX_SCREEN_OUTPUT_UPDATE_SIZE) {
				s = s.substring(0, MAX_SCREEN_OUTPUT_UPDATE_SIZE) + "\n" +
						CLIP_CLIP_LINE +
						"the rest of the screen output will be available when the job has finished";
				outputUpdateBufferFull = true;
			}
			updateScreenOutputCallback.accept(s);
		});

		// read process output stream
		BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {

				// check for script successful string
				if (successStringPattern.matcher(line).matches()) {
					if (!outputBufferFull) {
						screenOutput.append(line + "\n"); // include success string in the screen output for now
					}
					finishCallback.accept(JobState.RUNNING, screenOutput.toString());
					return;
				}

				// read normal output
				else {
					if (outputBufferFull) {
						continue;
					}

					//					// exculde print success string command
					//					if (line.contains(CompJob.SCRIPT_SUCCESSFUL_STRING)) {
					//						continue;
					//					}

					// make sure it always ends with \n
					line = line + "\n";

					// enough space in the buffer
					if (screenOutput.length() + line.length() <= MAX_SCREEN_OUTPUT_SIZE) {
						screenOutput.append(line);
						screenOutputSubject.onNext(true);
					} else {
						screenOutput.append(CLIP_CLIP_LINE);
						screenOutput.append("screen output was more than " + 
								MAX_SCREEN_OUTPUT_SIZE  + " characters long, the rest was discarded\n");

						outputBufferFull = true;
						screenOutputSubject.onNext(true);
					}
				}
			}
			
			// null line means end of stream --> job failed
			finishCallback.accept(JobState.FAILED, screenOutput.toString());
			return;

		} catch (IOException e) {
			// also canceling the job leads here 
			finishCallback.accept(JobState.ERROR, screenOutput.toString());
			return;
		} finally {
			screenOutputSubject.onComplete();
		}
	}

	public String getOutput() {
		return screenOutput.toString();
	}
}
