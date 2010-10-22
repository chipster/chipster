package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class Cigar {
	private List<CigarItem> elements = new ArrayList<CigarItem>();
	
	public void addElement(CigarItem e) {
		elements.add(e);
	}
	
	@Deprecated
	public long getReferenceIndex(long seqIndex) {
		
		long seqCounter = 0;
		long refCounter = 0;
		
		for (CigarItem element : elements) {
									
			if (element.getType().equals("M")) {
				
				seqCounter += element.getLength();
				refCounter += element.getLength();
				
			} else if (element.getType().equals("S")) {
				
				seqCounter += element.getLength();
				
			} else if (element.getType().equals("D") || 
					element.getType().equals("N") || 
					element.getType().equals("H")) {
				
				refCounter += element.getLength();
				
			} else if (element.getType().equals("I")) {
				
				seqCounter += element.getLength();
				
			} else if (element.getType().equals("P")) {
				
				//Do nothing, padded element exist only in some other read, not in this sequence nor in reference
			}
			
			if (seqCounter > seqIndex) {
				
				if (element.getType().equals("I")) {
					//There isn't reference sequence for insertion
					return -1; 
				} else if (element.getType().equals("S")) {
					//TODO how to draw this?
					return -1; 
				}
				
				return seqIndex + refCounter - seqCounter;
			}
		}
		
		//Request out of this read
		return -1;
	}
	
	/* 
	 * RegionContents are compared currently with the toString method, here we return anything
	 * constant so that comparison doesn't care about this object.
	 */
	@Override
	public String toString() {
		return "Cigar";
	}

	public static List<ReadPart> getVisibleRegions(RegionContent read) {
		Cigar cigar = (Cigar) read.values.get(ColumnType.CIGAR); // Cigar can be null

		if (cigar == null) {
			return Arrays.asList(new ReadPart[] { new ReadPart(read) });
			
		} else {
			LinkedList<ReadPart> regions = new LinkedList<ReadPart>();
			
			// Split read into regions using cigar
			long refCoord = read.region.start.bp;;
			long seqCoord = 0;
			String seq = (String) read.values.get(ColumnType.SEQUENCE);

			for (CigarItem element : cigar.elements) {
				
				if (element.isVisible()) {
					String subSeq = seq.substring((int)seqCoord, (int)(seqCoord + element.getLength()));
					ReadPart region = new ReadPart(refCoord, refCoord + element.getLength(), read.region.start.chr, subSeq);
					regions.add(region);
				}
				
				if (element.consumesReferenceBases()) {
					refCoord += element.getLength();
				}
				
				if (element.consumesReadBases()) {
					seqCoord += element.getLength();
				}

			}
			
			return regions;
		}
		
	}
}
