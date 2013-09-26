package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.CoverageTool;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;

public class CoverageStorage {
	
	private static final long AVERAGE_LIMIT = 2*1000;
	
	private TreeMap<BpCoord, Base> forwardBases = new TreeMap<BpCoord, Base>();
	private TreeMap<BpCoord, Base> reverseBases = new TreeMap<BpCoord, Base>();
	
	private TreeMap<Region, Float> forwardAverages = new TreeMap<Region, Float>();
	private TreeMap<Region, Float> reverseAverages = new TreeMap<Region, Float>();
	
	private boolean totalBasesNeedsRefresh;
	private boolean totalAveragesNeedsRefresh;
	
	private TreeMap<BpCoord, Base> totalBases;
	private TreeMap<Region, Float> totalAverages;

	
	public void addAverages(DataResult dataResult, Region filterRegion) {
		
		DataRequest request = dataResult.getRequest();

		if (request != null) {
			if (request.getRequestedContents().contains(DataType.COVERAGE_AVERAGE)) {

				CoverageTool.convertRegionContentListToFloatList(dataResult.getFeatures(), forwardAverages, reverseAverages);		
				filterAverages(forwardAverages, filterRegion);
				filterAverages(reverseAverages, filterRegion);
				totalAveragesNeedsRefresh = true;				
			}			
		}
	}

	public void addBaseCoverage(DataResult dataResult, Region filterRegion) {
				
		DataRequest request = dataResult.getRequest();

		if (request != null) {
			if (request.getRequestedContents().contains(DataType.COVERAGE)) {

				CoverageTool.convertRegionContentListToBaseList(dataResult.getFeatures(), forwardBases, reverseBases);				
				filterBases(forwardBases, filterRegion);
				filterBases(reverseBases, filterRegion);
				totalBasesNeedsRefresh = true;
			}
		}
	}

	private void filterBases(TreeMap<BpCoord, Base> bases,
			Region filterRegion) {
		
		Iterator<BpCoord> iterator = bases.keySet().iterator();
		
		while (iterator.hasNext()) {
			BpCoord region = iterator.next();
			
			if (!filterRegion.contains(region)) {
				iterator.remove();
			}
		}
	}

	private void filterAverages(TreeMap<Region, Float> averages, Region filterRegion) {

		Iterator<Region> iterator = averages.keySet().iterator();
		
		while (iterator.hasNext()) {
			Region region = iterator.next();
			
			if (!filterRegion.intersects(region)) {
				iterator.remove();
			}
		}
	}

	public TreeMap<BpCoord, Base> getTotalBases() {
		
		if (totalBases == null || totalBasesNeedsRefresh) {

			Iterator<Entry<BpCoord, Base>> forward = forwardBases.entrySet().iterator();		
			Iterator<Entry<BpCoord, Base>> reverse = reverseBases.entrySet().iterator();						 				

			totalBases = CoverageTool.getTotalBases(forward, reverse);
			totalBasesNeedsRefresh = false;
		}
		
		return totalBases;
	}
	
	public TreeMap<Region, Float> getTotalAverageCoverage() {
		
		if (totalAverages == null || totalAveragesNeedsRefresh) {
			
			Iterator<Entry<Region, Float>> forward = forwardAverages.entrySet().iterator();		
			Iterator<Entry<Region, Float>> reverse = reverseAverages.entrySet().iterator();
							
			totalAverages = CoverageTool.getTotalAverageCoverage(forward, reverse);
			
			totalAveragesNeedsRefresh = false;
		}
		
		return totalAverages;
	}


	public Base getBase(BpCoord location, Strand strand) {
		
		if (totalBases == null || totalBasesNeedsRefresh) {
			throw new IllegalStateException("totalBases is not initialized");
		}
		
		if (strand == Strand.FORWARD) {
			return forwardBases.get(location);
				
		} else if (strand == Strand.REVERSE) {				
			return reverseBases.get(location);
			
		} else if (strand == Strand.BOTH) {				
			return totalBases.get(location);
		}
		return null;
	}


	public Float getAverage(Region region, Strand strand) {
		
		if (totalAverages == null || totalAveragesNeedsRefresh) {
			throw new IllegalStateException("totalAverages is not initialized");
		}

		if (strand == Strand.FORWARD) {
			return forwardAverages.get(region);
				
		} else if (strand == Strand.REVERSE) {				
			return reverseAverages.get(region);
			
		} else if (strand == Strand.BOTH) {
			return totalAverages.get(region);
		}
		return null;
	}


	public boolean isAverage(Region viewRegion) {
		return viewRegion.getLength() >= AVERAGE_LIMIT;
	}

}
