package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.CigarElement;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CigarItem.CigarItemType;

/**
 * Represents the CIGAR string, as supported by the SAM standard. CIGAR string is used to
 * store the information on what parts of the read matched the genome. 
 * 
 * @author Aleksi Kallio, Petri Klemelä
 *
 */
public class Cigar {
	private List<CigarItem> cigarItems = new ArrayList<CigarItem>();
	private LinkedList<ReadPart> visibleElements = null;
	private LinkedList<ReadPart> processedElements = null;
	private Feature read;
	private net.sf.samtools.Cigar samCigar;

	public Cigar(Feature read, List<CigarItem> cigarItems) {
		this.read = read;
		this.cigarItems = cigarItems;
	}

	public Cigar(Feature read, net.sf.samtools.Cigar samCigar) {
		this.read = read;
		this.samCigar = samCigar;

		for (CigarElement picardElement : samCigar.getCigarElements()) {
			cigarItems.add(new CigarItem(picardElement));
		}
	}

	public void addElement(CigarItem e) {
	}

	@Deprecated
	public long getReferenceIndex(long seqIndex) {

		long seqCounter = 0;
		long refCounter = 0;

		for (CigarItem element : cigarItems) {

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
	public static List<ReadPart> splitElements(Feature read) {
		Cigar cigar = (Cigar) read.values.get(DataType.CIGAR); // Cigar can be null

		if (cigar == null) {
			return Arrays.asList(new ReadPart[] { new ReadPart(read) });

		} else {

			return cigar.getElements();
		}
	}

	/**
	 * Method for splitting a read according to defined CigarItemTypes. Sequence, cigar and region
	 * are splitted and only relevant parts are stored into result reads. Collection of all other
	 * values is cloned for each result read. Results don't include the splitters and corresponding 
	 * parts of the sequence and region. Returns the original read if there is no cigar. 
	 * 
	 * @param read
	 * @param splitters
	 * @return
	 */
	public static List<Feature> splitRead(Feature read, Collection<CigarItemType> splitters) {
		Cigar cigar = (Cigar) read.values.get(DataType.CIGAR); // Cigar can be null

		if (cigar == null) {
			return Arrays.asList(new Feature[] { read });

		} else {

			return cigar.splitReadWithCigar(read, splitters);
		}
	}

	private List<Feature> splitReadWithCigar(Feature read, Collection<CigarItemType> splitters) {

		List<Feature> splittedReads = new LinkedList<Feature>();

		long refCoord = read.region.start.bp;
		long seqCoord = 0;
		String seq = (String) read.values.get(DataType.SEQUENCE);

		String combinedSeq = "";
		List<CigarItem> combinedCigar = new LinkedList<CigarItem>();
		Region region = new Region(-1l, -1l, read.region.start.chr);


		for (CigarItem cigarItem : cigarItems) {
			String subSeq = "";
			if (cigarItem.consumesReadBases()) {
				subSeq = seq.substring((int)seqCoord, (int)(seqCoord + cigarItem.getLength()));
			}

			if (!splitters.contains(cigarItem.getCigarItemType())) {
				combinedSeq += subSeq;
				combinedCigar.add(cigarItem);
				if (region.start.bp == -1) {
					region.start.bp = refCoord;
				}
				region.end.bp = refCoord + cigarItem.getLength();
			} else {
				LinkedHashMap<DataType, Object> values = (LinkedHashMap<DataType, Object>) read.values.clone();
				Feature splittedRead = new Feature(region, values);
				values.put(DataType.CIGAR, new Cigar(splittedRead, combinedCigar));
				values.put(DataType.SEQUENCE, combinedSeq);
				splittedReads.add(splittedRead);

				combinedCigar = new LinkedList<CigarItem>();
				combinedSeq = "";
				region = new Region(-1l, -1l, read.region.start.chr);
			}

			if (cigarItem.consumesReferenceBases()) {
				refCoord += cigarItem.getLength();
			}

			if (cigarItem.consumesReadBases()) {
				seqCoord += cigarItem.getLength();
			}
		}

		LinkedHashMap<DataType, Object> values = (LinkedHashMap<DataType, Object>) read.values.clone();
		Feature splittedRead = new Feature(region, values);
		values.put(DataType.CIGAR, new Cigar(splittedRead, combinedCigar));
		values.put(DataType.SEQUENCE, combinedSeq);
		splittedReads.add(splittedRead);

		return splittedReads;
	}

	private List<ReadPart> getElements() {

		// Regions are parsed lazily
		if (processedElements == null) {
			processedElements = new LinkedList<ReadPart>();

			// Split read into regions using cigar
			long refCoord = read.region.start.bp;
			long seqCoord = 0;
			String seq = (String) read.values.get(DataType.SEQUENCE);

			for (CigarItem cigarItem : cigarItems) {
				String subSeq = null;
				if (cigarItem.consumesReadBases()) {
					subSeq = seq.substring((int)seqCoord, (int)(seqCoord + cigarItem.getLength()));
				}
				ReadPart region = new ReadPart(refCoord, refCoord + cigarItem.getLength(), read.region.start.chr, read, subSeq);
				region.setCigarItem(cigarItem);
				processedElements.add(region);

				if (cigarItem.consumesReferenceBases()) {
					refCoord += cigarItem.getLength();
				}

				if (cigarItem.consumesReadBases()) {
					seqCoord += cigarItem.getLength();
				}

			}
		}
		return processedElements;
	}

	public String toInfoString() {
		String str = "";

		if (samCigar != null) {
			str = samCigar.toString();
		} else {

			for (CigarItem item : cigarItems) {
				str += item.getLength() + item.getType();
			}
		}

		return "Cigar: " + str;
	}
}
