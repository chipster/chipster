package fi.csc.chipster.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.gbrowser.intervals.IntervalOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.IOUtils;

public class IntervalTool extends JavaAnalysisJobBase {

	@Override
	public String getSADL() {
		return 	"TOOL \"Interval operations\" / IntervalTool.java: \"Remove intervals\" (Removes intervals of the second dataset from the first dataset.)" + "\n" +
				"INPUT data1.bed: \"First input file\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Second input file\" TYPE GENERIC" + "\n" +
				"OUTPUT result.bed: \"Result file\"" + "\n" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs\" INTEGER FROM 1 DEFAULT 1 (How many basepairs are required for overlapping)";
	}
	
	
	@Override
	protected void execute() { 
		try {
			updateStateToClient(JobState.RUNNING, "preprocessing");
			File inputFile1 = new File(jobWorkDir, analysis.getInputFiles().get(0).getFileName()); 
			File inputFile2 = new File(jobWorkDir, analysis.getInputFiles().get(1).getFileName());
			
			IntervalOperations tool = new IntervalOperations();

			List<RegionContent> rows1 = tool.loadFile(inputFile1);
			List<RegionContent> rows2 = tool.loadFile(inputFile2);
			Long minOverlap = Long.parseLong(inputMessage.getParameters().get(0));

			LinkedList<BpCoordRegion> intersections = tool.intersect(rows1, rows2, minOverlap, false);

			FileOutputStream outputStream = null;
			try {
				outputStream = new FileOutputStream(new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID())); 
				tool.print(intersections, outputStream);

			} finally {
				IOUtils.closeIfPossible(outputStream);
			}
			
		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
			return;
		}
		updateStateToClient(JobState.RUNNING, "preprocessing finished");
	}
	


}
