package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Graphics;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;

import javax.swing.JFrame;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Cytoband;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.CytobandHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class CytobandDemo extends JFrame implements AreaResultListener {

	private static final int THICKNESS = 11;
	private static final int MARGIN = 50;
	
	private static final int WIN_WIDTH = 1280;

	CytobandHandlerThread dataThread;
	public Queue<AreaRequest> areaRequestQueue = new ConcurrentLinkedQueue<AreaRequest>();

	private SortedSet<Cytoband> cbands = new TreeSet<Cytoband>();

	final Map<Cytoband.Stain, Color> stainColors = new HashMap<Cytoband.Stain, Color>();

	{
		stainColors.put(Cytoband.Stain.GNEG, Color.white);
		stainColors.put(Cytoband.Stain.GPOS25, Color.lightGray);
		stainColors.put(Cytoband.Stain.GPOS33, Color.lightGray);
		stainColors.put(Cytoband.Stain.GPOS50, Color.gray);
		stainColors.put(Cytoband.Stain.GPOS66, Color.darkGray);
		stainColors.put(Cytoband.Stain.GPOS75, Color.darkGray);
		stainColors.put(Cytoband.Stain.GPOS100, Color.black);
		stainColors.put(Cytoband.Stain.GPOS, Color.black);
		stainColors.put(Cytoband.Stain.ACEN, null);
		stainColors.put(Cytoband.Stain.GVAR, Color.lightGray);
		stainColors.put(Cytoband.Stain.STALK, null);
		stainColors.put(Cytoband.Stain.TIP, null);
		stainColors.put(Cytoband.Stain.UNRECOGNIZED, null);
	}

	public CytobandDemo() {

		//Init Chipster data layer
		CytobandDataSource file = null;
		try {
			
			/*
			 *  Download and extract following files
			 *  ftp://ftp.ensembl.org/pub/release-65/mysql/rattus_norvegicus_core_65_34/seq_region.txt.gz
			 *  ftp://ftp.ensembl.org/pub/release-65/mysql/rattus_norvegicus_core_65_34/karyotype.txt.gz
			 *  
			 *  and adjust this paths correspondingly:
			 */
			
			file = new CytobandDataSource(new File("/home/klemela/chipster/ohtu/rattus/karyotype.txt"), 
					new File("/home/klemela/chipster/ohtu/rattus/seq_region.txt"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		dataThread = new CytobandHandlerThread(file, areaRequestQueue, this);
		dataThread.start();
		
		requestData();
	}
	
	public void requestData() {
		
		areaRequestQueue.add(new AreaRequest(
				new Region(0l, 270000000l, new Chromosome("1")), 
				new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {ColumnType.VALUE })),
				new FsfStatus()));
		
		areaRequestQueue.add(new AreaRequest(
				new Region(0l, 270000000l, new Chromosome("2")), 
				new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {ColumnType.VALUE })),
				new FsfStatus()));
		
		areaRequestQueue.add(new AreaRequest(
				new Region(0l, 270000000l, new Chromosome("X")), 
				new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {ColumnType.VALUE })),
				new FsfStatus()));

		dataThread.notifyAreaRequestHandler();
	}

	public int bpToDisplay(long bp) {
		return (int) (bp / 270000000f * WIN_WIDTH);
	}


	public void paint(Graphics g) {
		super.paint(g);

		//To determine if Stain.ACEN should be drawn "<" or ">"
		boolean firstGap = true;

		Iterator<Cytoband> cbandIter = cbands.iterator();
		
		//Only used to put different chromosomes to different rows
		Set<Chromosome> chrs = new HashSet<Chromosome>();
		
		while (cbandIter.hasNext()) {

			Cytoband cband = cbandIter.next();
			
			Cytoband.Stain stain = cband.getStain();
			Color stainColor = stainColors.get(stain);
			String text = (String) cband.getBand();
			
			chrs.add(cband.getRegion().start.chr);
			
			int x = bpToDisplay(cband.getRegion().start.bp);
			int width = bpToDisplay(cband.getRegion().getLength());	
			int y = chrs.size() * (MARGIN + THICKNESS * 3) + MARGIN;

			//Stain color is defined only for those types that are visualised with rectangle
			if (stainColor != null) {

				g.setColor(stainColor);			
				g.fillRect(x, y, width, THICKNESS);

				g.setColor(Color.black);		
				g.drawRect(x, y, width, THICKNESS);
				
				
				final int CHAR_WIDTH = 7;
				
				if (width > text.length() * CHAR_WIDTH) {

					g.setColor(Color.black);
					g.drawString(text, x, y + THICKNESS * 2);
				}

				firstGap = true;

			} else if (stain == Cytoband.Stain.ACEN) {

				int sideX = bpToDisplay(cband.getRegion().end.bp);
				int cornerX = bpToDisplay(cband.getRegion().start.bp);

				if (firstGap) {
					int tmp = sideX;
					sideX = cornerX;
					cornerX = tmp;
					firstGap = false;
				}

				g.setColor(Color.black);				
				g.drawLine(sideX, y, cornerX, y + THICKNESS / 2);
				g.drawLine(sideX, y + THICKNESS - 1, cornerX, y + THICKNESS / 2);

			} 
		}
	}

	public static void main(String arg[]) {

		CytobandDemo frame = new CytobandDemo();

		frame.setSize(WIN_WIDTH, 480);

		frame.setVisible(true);
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {
				
		for (RegionContent content : areaResult.getContents()) {

			Cytoband cband = (Cytoband) content.values.get(ColumnType.VALUE);
			cbands.add(cband);
		}
		
		this.repaint();
	}
}

