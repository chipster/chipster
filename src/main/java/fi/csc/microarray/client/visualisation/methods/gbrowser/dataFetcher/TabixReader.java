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

		long[][] lengthToLevel = new long[7][2];
		lengthToLevel[0][0] = 0L; lengthToLevel[0][1] = 2L;
		lengthToLevel[1][0] = 5000L; lengthToLevel[1][1] = 4L;
		lengthToLevel[2][0] = 10000L; lengthToLevel[2][1] = 8L;
		lengthToLevel[3][0] = 30000L; lengthToLevel[3][1] = 16L;
		lengthToLevel[4][0] = 50000L; lengthToLevel[4][1] = 32L;
		lengthToLevel[5][0] = 100000L; lengthToLevel[5][1] = 64L;
		lengthToLevel[6][0] = 1000000L; lengthToLevel[6][1] = 128L;

		String s = null;
		List<RegionContent> results = new LinkedList<RegionContent>();
		
		long length = Long.parseLong(endBp) - Long.parseLong(startBp);
		
		// hack the generic part of name
		String prefix = file;
		while (Character.isDigit(prefix.charAt(prefix.length()-1))) {
			prefix = prefix.substring(0, prefix.length()-1);
		}

		// decide summary level 
		long level = lengthToLevel[0][1];
		for (int i = 1; i < lengthToLevel.length; i++) {
			if (length > lengthToLevel[i][0]) {
				level = lengthToLevel[i][1];
			}
		}
		file = prefix + "-summary" + level;
		
		//FIXME replace with real Java library for tabix, for example 
		//https://genomeview.svn.sourceforge.net/svnroot/genomeview/jannot/trunk/net/sf/jannot/tabix/
		Process p = Runtime.getRuntime().exec(
				"/home/" + System.getProperty("user.name") + "/tabix-0.2.2/tabix " + 
				file + " " + 
				chr + ":" + 
				startBp + "-" + endBp);

		BufferedReader stdInput = new BufferedReader(new 
				InputStreamReader(p.getInputStream()));

		BufferedReader stdError = new BufferedReader(new 
				InputStreamReader(p.getErrorStream()));

		// read the output from the command
		while ((s = stdInput.readLine()) != null) {
			
			String[] splitted = s.split("\t");

			//Chromosome chr = splitted[0];
			Chromosome lineChr = new Chromosome("6");
			long start = Long.parseLong(splitted[1]);
			long end = Long.parseLong(splitted[2]);

			float value = Float.parseFloat(splitted[3]);

			Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
			values.put(ColumnType.VALUE, new Float(value));

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
