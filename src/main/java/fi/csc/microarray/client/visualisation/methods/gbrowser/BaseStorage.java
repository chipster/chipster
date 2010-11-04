package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class BaseStorage {

	private static final int MIN_SIGNIFICANT_SNP_COUNT = 2;
	private static final double MIN_SIGNIFICANT_SNP_RATIO = 0.20;

	private Object baseCache = null;
	private TreeMap<Long, Base> collector; 
	
	public enum Acid { 
		A, 
		C, 
		G, 
		T;

		public Acid complement() {
			switch (this) {
				case A:
					return T;
				case C:
					return G;
				case G:
					return C;
				case T:
					return A;
				default:
					throw new RuntimeException("this never happens");
			}
		}

		public static Acid fromCharacter(char character) {
			switch (Character.toUpperCase(character)) {
			case 'A':
				return A;
			case 'C':
				return C;
			case 'G':
				return G;
			case 'T':
				return T;
			default:
				return null;
			}
		}

		public static Acid fromCharacter(char character, boolean complement) {
			Acid acid = fromCharacter(character);
			
			if (acid != null && complement) {
				return acid.complement();
			} else {
				return acid;
			}
		} 
	};


	public class Base {

		private int[] acidCounts = new int[Acid.values().length];
		private Long bpLocation;
		private int[] snpCounts = null;
		private int totalSNPCount = 0;
		private int totalCount = 0;
		private Acid referenceAcid;
		
		public Base(Long bpLocation, Acid referenceAcid) {
			
			this.bpLocation = bpLocation;
			this.referenceAcid = referenceAcid;
		}

		public Long getBpLocation() {
			return bpLocation;
		}
		
		public int getCoverage() {
			return totalCount;
		}
		
		public boolean hasSignificantSNPs() {
			int[] snpCounts = getSNPCounts();

			// Try to find good enough SNP acid
			for (Acid acid : Acid.values()) {
				if (snpCounts[acid.ordinal()] > 0) {
					return true;
				}
//				if (snpCounts[acid.ordinal()] >= MIN_SIGNIFICANT_SNP_COUNT) {
//					if (((double)snpCounts[acid.ordinal()])/((double)totalSNPCount) >= MIN_SIGNIFICANT_SNP_RATIO) {
//						return true;
//					}
//				}
			}
			
			// None found
			return false;
		}

		public int[] getAcidCounts() {
			return acidCounts;
		}

		public void addAcid(Acid acid) {
			if (snpCounts != null) {
				throw new IllegalStateException("cannot add acids after SNP counts have been calculated");
			}
			acidCounts[acid.ordinal()]++;
			totalCount++;
		}
		
		public int getTotalSNPCount() {
			if (snpCounts == null) {
				throw new IllegalStateException("cannot get total before calling getSNPCounts()");
			}
			return totalSNPCount;
		}
		
		public int[] getSNPCounts() {
			if (snpCounts == null) {
				snpCounts = new int[Acid.values().length];
				
				// Mark SNP's, if possible
				if (referenceAcid != null) {
					for (Acid acid : Acid.values()) {
						if (acid.compareTo(referenceAcid) == 0) {
							snpCounts[acid.ordinal()] = 0;
						} else {
							snpCounts[acid.ordinal()] = acidCounts[acid.ordinal()];
							totalSNPCount += acidCounts[acid.ordinal()];
						}
					}
				}
			}
			return snpCounts;
		}
	}

	public void updateBases(Collection<RegionContent> reads, View view) {

		
		// Divide to two groups: in view and out of view
		LinkedList<RegionContent> readsInView = new LinkedList<RegionContent>();
		LinkedList<RegionContent> readOutOfView = new LinkedList<RegionContent>();
		Iterator<RegionContent> iter = reads.iterator();
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			if (read.region.intercepts(view.getBpRegion())) {
				readsInView.add(read);
			} else {
				readOutOfView.add(read);
			}
		}
		
		// Update reads
		// TODO
		
		// Decide base cache update strategy and execute it
		if (baseCache != null && readsInView.size() > readOutOfView.size()) {
			
			// Cache exists and most of the reads still in, so update old
			for (RegionContent toBeRemoved : readOutOfView) {
				removeRead(toBeRemoved);
			}
			
		} else {
			
			// Most of the reads are out, so regenerate everything
			baseCache = new Object();
			for (RegionContent read : readsInView) {
				addRead(read);
			}
		}
		
		
	}

	private void addRead(RegionContent read) {
		// TODO Auto-generated method stub
	}

	private void removeRead(RegionContent toBeRemoved) {
		// TODO Auto-generated method stub
	}

	/**
	 * Goes through data and gives count for each location and acid.
	 * @param refSeq 
	 */
	public void getAcidCounts(Collection<RegionContent> reads, View view, char[] refSeq) {
	
		// Sweep collector
		collector = new TreeMap<Long, Base>();
		Iterator<RegionContent> iter = reads.iterator();

		// iterate over RegionContent objects (one object corresponds to one read)
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// remove those that are not in this view
			if (!read.region.intercepts(view.getBpRegion())) {
				iter.remove();
				continue;
			}

			Strand strand = (Strand) read.values.get(ColumnType.STRAND);

			for (ReadPart readPart : Cigar.splitVisibleElements(read)) {

				// Skip elements that are not in this view
				if (!readPart.intercepts(view.getBpRegion())) {
					continue;
				}

				Base base = null;

				String seq = readPart.getSequencePart();
				for (int j = 0; j < seq.length(); j++) {

					Long bp = readPart.start.bp + j;
					
					// Part of read can be out of view
					if (bp.longValue() < view.bpRegion.start.bp.longValue()) {
						continue;
						
					} else if (bp.longValue() > view.bpRegion.end.bp.longValue()) {
						break;
					}

					if (!collector.containsKey(bp)) {
						
						int viewIndex = bp.intValue() - view.bpRegion.start.bp.intValue();
						Acid referenceAcid = Acid.fromCharacter(refSeq[viewIndex], strand == Strand.REVERSED);
						
						base = new Base(bp, referenceAcid);
						collector.put(bp, base);
						
					} else {
						base = collector.get(bp);
					}

					Acid acid = Acid.fromCharacter(seq.charAt(j), strand == Strand.REVERSED);
					if (acid != null) {
						base.addAcid(acid);
					}
				}
			}
		}
	}

	public Iterator<Base> iterator() {
		return collector.values().iterator();
	}
}
