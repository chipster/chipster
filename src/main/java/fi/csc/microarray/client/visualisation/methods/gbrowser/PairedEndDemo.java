package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Graphics;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Queue;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;

import javax.swing.JFrame;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.SAMHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class PairedEndDemo extends JFrame implements AreaResultListener {

	private static final int WIN_WIDTH = 1280;

	public SAMHandlerThread dataThread;
	public Queue<AreaRequest> areaRequestQueue = new ConcurrentLinkedQueue<AreaRequest>();

	private SortedSet<RegionContent> reads = new TreeSet<RegionContent>();

	public PairedEndDemo() {

		//Init Chipster data layer
		SAMDataSource file = null;
		try {
			
			// Adjust these paths to point to the demo data
			
			file = new SAMDataSource(new File("ohtu-between-chrs.bam"),
					new File("ohtu-between-chrs.bam.bai"));
			
//			file = new SAMDataSource(new File("ohtu-paired-end.bam"),
//					new File("ohtu-paired-end.bam.bai"));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		dataThread = new SAMHandlerThread(file, areaRequestQueue, this);
		dataThread.start();
		
		requestData();
	}
	
	public void requestData() {
		
		/* This is the best place for limiting the amount of data. If all connections make this too slow,
		 * request only data of smaller areas. 
		 */
		
		areaRequestQueue.add(new AreaRequest(
				new Region(0l, 270000000l, new Chromosome("1")), 
				new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
						ColumnType.ID, ColumnType.STRAND, ColumnType.MATE_POSITION})),
				new FsfStatus()));
		
		areaRequestQueue.add(new AreaRequest(
				new Region(0l, 270000000l, new Chromosome("X")), 
				new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
						ColumnType.ID, ColumnType.STRAND, ColumnType.MATE_POSITION})),
				new FsfStatus()));

		dataThread.notifyAreaRequestHandler();
	}

	public int bpToDisplay(long bp) {
		return (int) (bp / 270000000f * WIN_WIDTH);
	}

 
	public void paint(Graphics g) {
		
		//Background painting disabled to avoid flickering, demo needs to paint only once
		//super.paint(g);

		paint(g, new Chromosome("1"), 50, new Chromosome("X"), 400, Color.black, Color.green);
		paint(g, new Chromosome("X"), 400, new Chromosome("1"), 50, Color.black, Color.blue);
		
	}
	
	public void paint(Graphics g, Chromosome fromChr, int fromY, Chromosome toChr, int toY, Color readColor, Color connectionColor) {
		
		
		Iterator<RegionContent> readIter = reads.iterator();
		
		while (readIter.hasNext()) {

			RegionContent read = readIter.next();
			BpCoord mate = (BpCoord)read.values.get(ColumnType.MATE_POSITION);
			
			if (!fromChr.equals(read.region.start.chr) || !toChr.equals(mate.chr)) {
				continue;
			}
			
			int x = bpToDisplay(read.region.start.bp);
			// End coordinate of the read is read.region.end
			
			
			g.setColor(readColor);
			g.drawRect(x, fromY, 10, 10);
			
			int mateX = bpToDisplay(mate.bp);
			
			g.setColor(connectionColor);
			g.drawLine(mateX, toY, x, fromY);
		}
	}

	public static void main(String arg[]) {

		PairedEndDemo frame = new PairedEndDemo();

		frame.setSize(WIN_WIDTH, 480);

		frame.setVisible(true);
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {
		
		Chromosome readChr = null;
		Chromosome mateChr = null;
		
		Chromosome upperChr = new Chromosome("1");
		Chromosome lowerChr = new Chromosome("X");
		Chromosome unmappedChr = new Chromosome("*");
		
				
		for (RegionContent read : areaResult.getContents()) {
			
			readChr = read.region.start.chr;
			
			if (!upperChr.equals(readChr) && !lowerChr.equals(readChr)) {
				//Make sure that there aren't any data from other chromosomes
				continue;
			}
			
			mateChr = ((BpCoord)read.values.get(ColumnType.MATE_POSITION)).chr;
			
			if (unmappedChr.equals(mateChr)) {
				//Mate is unmapped, location is unknown, 
				//These are removed already from the demo data
				continue;
			}
			
			if (readChr.equals(mateChr)) {
				
				//Mate is in the same chromosome, not visualised in this demo
				continue;
			}
			
			
			if (!upperChr.equals(mateChr) && !lowerChr.equals(mateChr)) {
				
				//Only showing connections between Chromosome 1 and X in this demo, reject others
				continue;
			}
			
			/*
			 * The connections of read pairs aren't stored in the data separately, but both ends 
			 * know the location of other pair: 
			 * 
			 * (BpCoord)read.values.get(ColumnType.MATE_POSITION)
			 * 
			 * If there are more than one read in that location, the correct one can be found by
			 * comparing the identifiers of the reads, pairs have identical identifiers: 
			 * 
			 *(String) read.values.get(ColumnType.ID)
			 * 
			 */
			
			reads.add(read);
			
		}
		
		this.repaint();
	}
}

