package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;

public class CoverageTool {
	
	public static final int BIN_SIZE = 16;

	public static TreeMap<BpCoord, Base> getTotalBases(Iterator<Entry<BpCoord, Base>> forward, Iterator<Entry<BpCoord, Base>> reverse) {
		TreeMap<BpCoord, Base> totalBases = new TreeMap<BpCoord, Base>();
				
		addBaseCounts(forward, totalBases);		
		addBaseCounts(reverse, totalBases);
				
		return totalBases;
	}
	
	public static TreeMap<Region, Float> getTotalAverageCoverage(
			Iterator<Entry<Region, Float>> forward,
			Iterator<Entry<Region, Float>> reverse) {
		
		TreeMap<Region, Float> totals = new TreeMap<Region, Float>();
		
		addFloatCounts(forward, totals);
		addFloatCounts(reverse, totals);
		
		return totals;
		
	}
	
	private static void addBaseCounts(Iterator<Entry<BpCoord, Base>> iterator,
			TreeMap<BpCoord, Base> totalBases) {
		
		while (iterator.hasNext()) {
						
			Entry<BpCoord, Base> entry = iterator.next();
			
			Base totalBase = totalBases.get(entry.getKey());
						
			if (totalBase == null) {
						
				totalBase = new Base(entry.getKey().bp, null);
				totalBases.put(entry.getKey(), totalBase);
			}
			
			
			int[] counts = entry.getValue().getNucleotideCounts();
			int[] totalCounts = totalBase.getNucleotideCounts();
			
			for (int i = 0; i < totalCounts.length; i++) {
				totalCounts[i] += counts[i];
			}
			
			totalBase.setNucleotideCounts(totalCounts);			
		}
	}

	private static void addFloatCounts(Iterator<Entry<Region, Float>> iterator,
			TreeMap<Region, Float> totals) {
		
		while (iterator.hasNext()) {
			
			Entry<Region, Float> entry = iterator.next();
			
			Region region = entry.getKey();			
						
			if (!totals.containsKey(region)) {
				totals.put(region, 0f);
			}
			
			float value = totals.get(region);
			
			totals.put(region, value + entry.getValue());
		}
	}

	public static void convertRegionContentListToBaseList(
			List<Feature> regionContentList, TreeMap<BpCoord, Base> forwardBases,
			TreeMap<BpCoord, Base> reverseBases) {
		
		for (Feature regCont : regionContentList) {
			
			Object value = regCont.values.get(DataType.VALUE);

			if (value != null && value instanceof Base) {
				Base base = (Base) value;

				if (regCont.values.get(DataType.STRAND) == Strand.FORWARD) {
					forwardBases.put(regCont.region.start, base);

				} else if (regCont.values.get(DataType.STRAND) == Strand.REVERSE) {
					reverseBases.put(regCont.region.start, base);
				}
			}			
		}		
	}
	

	public static void convertRegionContentListToFloatList(
			List<Feature> regionContentList,
			TreeMap<Region, Float> forwardAverages,
			TreeMap<Region, Float> reverseAverages) {
		
		for (Feature regCont : regionContentList) {
			
			Object value = regCont.values.get(DataType.COVERAGE_AVERAGE);

			if (value != null && value instanceof Float) {
				Float floatValue = (Float)value;

				if (regCont.values.get(DataType.STRAND) == Strand.FORWARD) {
					forwardAverages.put(regCont.region, floatValue);

				} else if (regCont.values.get(DataType.STRAND) == Strand.REVERSE) {
					reverseAverages.put(regCont.region, floatValue);
				}
			}			
		}
	}

	public static TreeMap<Long, LinkedList<Base>> binBases (
			TreeMap<BpCoord, Base> bases) {		

		TreeMap<Long, LinkedList<Base>> bins = new TreeMap<Long, LinkedList<Base>>();		

		for (Base base : bases.values()) {

			Long bin = getBin(base.getBpLocation());

			if (!bins.containsKey(bin)) {
				bins.put(bin, new LinkedList<Base>());
			}

			bins.get(bin).add(base);
		}

		return bins;
	}

	public static long getBin(long location) {
		return ((long)(location / BIN_SIZE)) * BIN_SIZE;
	}
				
	public static LinkedList<Feature> average(TreeMap<Long, LinkedList<Base>> bins, Chromosome chr, Strand strand) {
			
		LinkedList<Feature> resultList = new LinkedList<Feature>();
		
		for (Entry<Long, LinkedList<Base>> entry : bins.entrySet()) {
			
			float average = 0;
			
			for (Base base : entry.getValue()) {
				
				average += base.getCoverage();
			}
			
			average /= entry.getValue().size();
				
			Region region = new Region(entry.getKey(), entry.getKey() + BIN_SIZE, chr);
				
			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
			values.put(DataType.COVERAGE_AVERAGE, average);
			values.put(DataType.STRAND, strand);

			resultList.add(new Feature(region, values));		
		}
		return resultList;
	}

	public static LinkedList<Feature> average(
			LinkedList<Feature> resultList, Chromosome chr) {
		
		TreeMap<BpCoord, Base> forwardBases = new TreeMap<BpCoord, Base>();
		TreeMap<BpCoord, Base> reverseBases = new TreeMap<BpCoord, Base>();
		
		CoverageTool.convertRegionContentListToBaseList(resultList, forwardBases, reverseBases);
		
		TreeMap<Long, LinkedList<Base>> forwardBins = CoverageTool.binBases(forwardBases);
		TreeMap<Long, LinkedList<Base>> reverseBins = CoverageTool.binBases(reverseBases);
		
		LinkedList<Feature> averageCoverage = new LinkedList<Feature>();
		averageCoverage.addAll(CoverageTool.average(forwardBins, chr, Strand.FORWARD));
		averageCoverage.addAll(CoverageTool.average(reverseBins, chr, Strand.REVERSE));
		
		return averageCoverage;
	}
}
