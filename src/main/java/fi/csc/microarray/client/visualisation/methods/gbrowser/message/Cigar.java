package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.CigarElement;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * Represents the CIGAR string, as supported by the SAM standard. CIGAR string is used to
 * store the information on what parts of the read matched the genome. 
 * 
 * @author Aleksi Kallio
 *
 */
public class Cigar {
	private List<CigarItem> elements = new ArrayList<CigarItem>();
	private LinkedList<ReadPart> visibleElements = null;
	private RegionContent read;
	private net.sf.samtools.Cigar samCigar;
	
	public Cigar(RegionContent read, net.sf.samtools.Cigar samCigar) {
		this.read = read;
		this.samCigar = samCigar;
		
    	for (CigarElement picardElement : samCigar.getCigarElements()) {
    		elements.add(new CigarItem(picardElement));
    	}
	}
	
	public void addElement(CigarItem e) {
	}
	
	// FIXME remove this and use instead ReadpartDataProvider
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

	/**
	 * Utility method for splitting a read if it has cigar, otherwise returns the whole
	 * read wrapped into ReadPart.
	 * 
	 * @param read read to split
	 * 
	 * @return list of ReadPart objects, one for each of the Cigar elements
	 */
	public static List<ReadPart> splitVisibleElements(RegionContent read) {
		Cigar cigar = (Cigar) read.values.get(ColumnType.CIGAR); // Cigar can be null

		if (cigar == null) {
			return Arrays.asList(new ReadPart[] { new ReadPart(read) });

		} else {
			return cigar.getVisibleElements();
		}

	}
	
	private List<ReadPart> getVisibleElements() {
		
		// Regions are parsed lazily
		if (visibleElements == null) {
			visibleElements = new LinkedList<ReadPart>();

			// Split read into regions using cigar
			long refCoord = read.region.start.bp;
			long seqCoord = 0;
			String seq = (String) read.values.get(ColumnType.SEQUENCE);

			for (CigarItem element : elements) {

				if (element.isVisible()) {
					String subSeq = seq.substring((int)seqCoord, (int)(seqCoord + element.getLength()));
					ReadPart region = new ReadPart(refCoord, refCoord + element.getLength(), read.region.start.chr, read, subSeq);
					visibleElements.add(region);
				}

				if (element.consumesReferenceBases()) {
					refCoord += element.getLength();
				}

				if (element.consumesReadBases()) {
					seqCoord += element.getLength();
				}

			}
		}
		return visibleElements;
	}
	
	public String toInfoString() {
		return "Cigar: " + samCigar.toString();
	}
}
