package fi.csc.microarray.client.dataimport;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.channels.ClosedByInterruptException;
import java.nio.channels.ClosedChannelException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import sun.nio.ch.ChannelInputStream;

/**
 * Class to analyse file which will be imported. While analysing file 
 * the analyser tries to guess which is the best suitable delimeter to 
 * chop the data of the file. It also count lines on file.
 * 
 * @author Mikko Koski, Petri Klemelä
 *
 */
public class FileAnalyser {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(FileAnalyser.class);

	/**
	 * Input file
	 */
	private File file;
	
	/**
	 * Best suitable delimeter the analyser found
	 */
	private Delimiter delimiter;
	
	/**
	 * Lines on file
	 */
	private int linesOnFile;

	/**
	 * Is the file analysed
	 */
	private boolean isAnalysed;
	
	/**
	 * Delimiter guess interval start value. The actual interval is not linear. It is 
	 * doubled during each iteration. So this is only the start value.
	 */
	private final static int DELIMETER_GUESS_INTERVAL_START = 10;

	public FileAnalyser(File file) {
		this.file = file;
		this.isAnalysed = false;
	}

	
	/**
	 * Starts analysing the file
	 * 
	 * @throws IOException the exception is thrown if there occures some problems 
	 * while reading the file but also if the analysis is interrupted
	 */
	public void startAnalysing() throws IOException{				
		
		// Let's initialize the delimeter map
		Map<Delimiter, Integer> hitCountsOnFile = new HashMap<Delimiter, Integer>();
		for(Delimiter delim : Delimiter.values()){
			hitCountsOnFile.put(delim, 0);
		}
		
		// Line counter
		int lines = 0;
		
		// Variable for line
		String nextLine;
		
		BufferedReader buf = null;
		int delimInterval = FileAnalyser.DELIMETER_GUESS_INTERVAL_START;
		try {
			buf = new BufferedReader(new InputStreamReader(new ChannelInputStream((new RandomAccessFile(file, "r")).getChannel())));
			
			// Read through the file line by line
			while((nextLine = buf.readLine()) != null && !Thread.interrupted()){
				lines++;
				if(lines % delimInterval == 0){
					
					// Guess the best suitable delimeter for this line
					Delimiter bestSuitable = FileAnalyser.guessDelimeterForLine(nextLine);
					
					// Increase the hit count
					int delimeterCount = hitCountsOnFile.get(bestSuitable);
					delimeterCount++;
					hitCountsOnFile.put(bestSuitable, delimeterCount);
					delimInterval = delimInterval * 2;
				}
			}
		} catch (ClosedByInterruptException cbie){
			if(buf != null){
				buf.close();
			}
			throw cbie;
		} catch (ClosedChannelException che){
			if(buf != null){
				buf.close();
			}
			throw che;
		} finally {
			if(buf != null){
				buf.close();
			}
		}
		
		// Iterate through the delimiter map and search the delimimer 
		// which has most best suitable hits
		Delimiter bestSuitable = null;
		int mostHits = 0;
		for(Delimiter delim : hitCountsOnFile.keySet()){
			if(hitCountsOnFile.get(delim) >= mostHits){
				bestSuitable = delim;
				mostHits = hitCountsOnFile.get(delim);
			}
		}
		
		logger.debug("File analysing ended. Total line count is " + lines + " and the best suitable delimeter is " + bestSuitable.getName());
		
		// Save the results
		this.delimiter = bestSuitable;
		this.linesOnFile = lines;
		this.isAnalysed = true;
	}
	
	
	
	/**
	 * Guesses the best delimeter for this given line. The line is splitted 
	 * with every delimeter and the delimeter which gives the most splits 
	 * is chosen to the best delimeter 
	 * 
	 * @param line line to be splitted
	 * @return the best delimeter
	 */
	private static Delimiter guessDelimeterForLine(String line) {
		
		// Map which keeps count of hits per delimeter
		Map<Delimiter, Integer> hitCountsOnLine = new HashMap<Delimiter, Integer>();
		
		String[] splittedLine;
		// Iterates through all delimeters and splits the line with each one
		for(Delimiter delim : Delimiter.values()){
			splittedLine = line.split(delim.toString());
			
			// Store the hit count to the map
			hitCountsOnLine.put(delim, splittedLine.length);
		}
		
		Delimiter bestSuitable = null;
		int mostHits = 0;
		
		// Iterate through the hit count map and find out which delimeter got the most 
		// splitted line
		for(Delimiter delim : hitCountsOnLine.keySet()){
			// logger.debug("Hits for " + delim.getName() + ": " + hitCount.get(delim));
			if(hitCountsOnLine.get(delim) >= mostHits){
				bestSuitable = delim;
				mostHits = hitCountsOnLine.get(delim);
			}
		}
		
		logger.debug("The best suitable delimeter for this line is " + bestSuitable.getName());
		return bestSuitable;
	}

	public Delimiter getDelimiter() {
		if(isAnalysed){
			return delimiter;
		} else {
			throw new IllegalStateException("File is not analysed yet");
		}
	}

	public int getLinesOnFile() {
		if(isAnalysed){
			return linesOnFile;
		} else {
			throw new IllegalStateException("File is not analysed yet");
		}
	}
	
	/**
	 * Sets new input file
	 * @param inputFile
	 */
	public void setInputFile(File inputFile){
		if(!isSameFileAs(inputFile)){
			this.file = inputFile;
			this.delimiter = null;
			this.linesOnFile = -1;
			this.isAnalysed = false;
		} else {
			// The file is same file and it is already analysed. 
			// No need to do anything
		}
	}
	
	/**
	 * Returns if the file is already analysed
	 * 
	 * @param file another file
	 * @return is file already analysed
	 */
	public boolean isAnalysed(){
		return this.isAnalysed;
	}
	
	public boolean isSameFileAs(File file){
		return this.file.equals(file);
	}
}
