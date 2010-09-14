package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TabixReader {

	public static List<RegionContent> query(String file, String chr, String startBp, String endBp) throws IOException {

		String s = null;
		List<RegionContent> results = new LinkedList<RegionContent>();
		
		
		//FIXME replace with real Java library for tabix, for example 
		//https://genomeview.svn.sourceforge.net/svnroot/genomeview/jannot/trunk/net/sf/jannot/tabix/
		Process p = Runtime.getRuntime().exec(
				"/home/klemela/tabix-0.2.2/tabix " + 
				file + " " + 
				chr + ":" + 
				startBp + "-" + endBp);

		BufferedReader stdInput = new BufferedReader(new 
				InputStreamReader(p.getInputStream()));

		BufferedReader stdError = new BufferedReader(new 
				InputStreamReader(p.getErrorStream()));

		// read the output from the command
		while ((s = stdInput.readLine()) != null) {

			System.out.println(" * " + s);
			
			String[] splitted = s.split("\t");

			//Chromosome chr = splitted[0];
			Chromosome lineChr = new Chromosome("6");
			long start = Long.parseLong(splitted[1]);
			long end = Long.parseLong(splitted[2]);

			int value = Integer.parseInt(splitted[3]);

			Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
			values.put(ColumnType.VALUE, new Integer(value));

			RegionContent reg = new RegionContent(
					new BpCoordRegion(start, end, lineChr), values);
			results.add(reg);
		}

		// read any errors from the attempted command
		while ((s = stdError.readLine()) != null) {
			System.out.println("Error in TabixReader: ");
			System.out.println(s);
		}

		return results;
	}
}
