package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Graphics;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;

import javax.swing.JFrame;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class BedDemo extends JFrame implements AreaResultListener {

	private static final int THICKNESS = 11;
	private static final int MARGIN = 50;

	private static final int WIN_WIDTH = 1280;

	ChunkTreeHandlerThread dataThread;
	public Queue<AreaRequest> areaRequestQueue = new ConcurrentLinkedQueue<AreaRequest>();

	private SortedSet<Region> regions = new TreeSet<Region>();


	public BedDemo() {

		/* 
		 *  adjust paths:
		 */

		String dataPath = System.getProperty("user.home") + "/chipster/ohtu/";
		File BED_FILE = new File(dataPath + "peaks.bed");

		DataSource bedDataSource = null;
		try {
			bedDataSource = new ChunkDataSource(BED_FILE, new BEDParser());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		dataThread = new ChunkTreeHandlerThread(bedDataSource, areaRequestQueue, this);
		dataThread.start();

		requestData();
	}

	public void requestData() {

		areaRequestQueue.add(new AreaRequest(
				new Region(0l, Long.MAX_VALUE, new Chromosome("1")), 
				new HashSet<ColumnType>(),
				new FsfStatus()));

		areaRequestQueue.add(new AreaRequest(
				new Region(0l, Long.MAX_VALUE, new Chromosome("2")), 
				new HashSet<ColumnType>(),
				new FsfStatus()));

		areaRequestQueue.add(new AreaRequest(
				new Region(0l, Long.MAX_VALUE, new Chromosome("X")), 
				new HashSet<ColumnType>(),
				new FsfStatus()));

		dataThread.notifyAreaRequestHandler();
	}

	public int bpToDisplay(long bp) {
		return (int) (bp / 270000000f * WIN_WIDTH); //scale to show full chromosome, 270000000 is about the length of longest one
	}


	public void paint(Graphics g) {
		super.paint(g);

		Iterator<Region> regionIter = regions.iterator();

		//Only used to put different chromosomes to different rows
		Set<Chromosome> chrs = new HashSet<Chromosome>();

		while (regionIter.hasNext()) {

			Region region = regionIter.next();

			chrs.add(region.start.chr);

			int x = bpToDisplay(region.start.bp);
			int width = bpToDisplay(region.getLength());	
			int y = chrs.size() * (MARGIN + THICKNESS * 3) + MARGIN;

			g.setColor(Color.green.darker());			
			g.fillRect(x, y, width, THICKNESS);

			g.setColor(Color.black);		
			g.drawRect(x, y, width, THICKNESS);
		}
	}

	public static void main(String arg[]) {

		BedDemo frame = new BedDemo();

		frame.setSize(WIN_WIDTH, 480);

		frame.setVisible(true);
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent content : areaResult.getContents()) {

			regions.add(content.region);
		}

		this.repaint();
	}
}

