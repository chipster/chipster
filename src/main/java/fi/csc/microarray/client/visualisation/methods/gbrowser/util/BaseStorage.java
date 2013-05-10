package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Utility class for storing information that is organised by bases (genomic coordinates).
 * Helpful when you need to e.g. convert a set or reads to base specific coverage.
 * 
 * @author Aleksi Kallio
 *
 */
public class BaseStorage {

	private static final int MIN_SIGNIFICANT_SNP_COUNT = 2;
	private static final double MIN_SIGNIFICANT_SNP_RATIO = 0.20;

	private TreeMap<Long, Base> collector; 
	
	public enum Nucleotide { 
		A, 
		C, 
		G, 
		T;

		public Nucleotide complement() {
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

		public static Nucleotide fromCharacter(char character) {
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

		public static Nucleotide fromCharacter(char character, boolean complement) {
			Nucleotide nt = fromCharacter(character);
			
			if (nt != null && complement) {
				return nt.complement();
			} else {
				return nt;
			}
		} 
	};


	public static class Base {

		private int[] nucleutideCounts = new int[Nucleotide.values().length];
		private Long bpLocation;
		private int[] snpCounts = null;
		private int totalSNPCount = 0;
		private int totalCount = 0;
		private Nucleotide referenceNt;
		
		public Base(Long bpLocation, Nucleotide referenceNt) {
			this.bpLocation = bpLocation;
			this.referenceNt = referenceNt;
		}

		public Long getBpLocation() {
			return bpLocation;
		}
		
		public int getCoverage() {
			return totalCount;
		}
		
		public boolean hasSignificantSNPs() {
			int[] snpCounts = getSNPCounts();

			// Try to find good enough SNP nucleotide
			for (Nucleotide nt : Nucleotide.values()) {
				if (snpCounts[nt.ordinal()] >= MIN_SIGNIFICANT_SNP_COUNT) {
					if (((double)snpCounts[nt.ordinal()])/((double)totalCount) >= MIN_SIGNIFICANT_SNP_RATIO) {
						return true;
					}
				}
			}
			
			// None found
			return false;
		}

		public int[] getNucleotideCounts() {
			return nucleutideCounts;
		}

		public void addNucleotidee(Nucleotide nt) {
			if (snpCounts != null) {
				throw new IllegalStateException("cannot add nucleotides after SNP counts have been calculated");
			}
			nucleutideCounts[nt.ordinal()]++;
			totalCount++;
		}
		
		public int getTotalSNPCount() {
			if (snpCounts == null) {
				throw new IllegalStateException("cannot get total before calling getSNPCounts()");
			}
			return totalSNPCount;
		}
		
		public Nucleotide getReferenceNucleotide() {
			return referenceNt;
		}
		
		public int[] getSNPCounts() {
			if (snpCounts == null) {
				snpCounts = new int[Nucleotide.values().length];
				
				// Mark SNP's, if possible
				if (referenceNt != null) {
					for (Nucleotide nt : Nucleotide.values()) {
						if (nt.compareTo(referenceNt) == 0) {
							snpCounts[nt.ordinal()] = 0;
						} else {
							snpCounts[nt.ordinal()] = nucleutideCounts[nt.ordinal()];
							totalSNPCount += nucleutideCounts[nt.ordinal()];
						}
					}
				}
			}
			return snpCounts;
		}

		public void setNucleotideCounts(int[] nucleotideCounts) {
			this.nucleutideCounts = nucleotideCounts;
			
			totalCount = 0;
			for (int i = 0; i < nucleotideCounts.length; i++) {
				totalCount += nucleotideCounts[i];
			}
		}
	}

	/**
	 * Goes through data and gives count for each location and nucleotide.
	 * @param refSeq 
	 */
	public void getNucleotideCounts(Iterable<ReadPart> readParts, GBrowserView view, char[] refSeq) {
	
		// Sweep collector
		collector = new TreeMap<Long, Base>();

			for (ReadPart readPart : readParts) {

				// Skip elements that are not in this view
				if (view != null && !view.requestIntersects(readPart)) {
					continue;
				}
				
				// Skip invisible types
				if (!readPart.isVisible()) {
					continue;
				}

				Base base = null;

				String seq = readPart.getSequencePart();
				for (int j = 0; j < seq.length(); j++) {

					Long bp = readPart.start.bp + j;
					
					// Part of read can be out of view
					if (view != null && bp.longValue() < view.bpRegion.start.bp.longValue()) {
						continue;
						
					} else if (view != null && bp.longValue() > view.bpRegion.end.bp.longValue()) {
						break;
					}

					if (!collector.containsKey(bp)) {
						
						if (view != null && refSeq != null) {
							int viewIndex = bp.intValue() - view.bpRegion.start.bp.intValue();
							Nucleotide referenceNucleotide = Nucleotide.fromCharacter(refSeq[viewIndex]);						
							base = new Base(bp, referenceNucleotide);
						} else {
							base = new Base(bp, null);
						}
						
						collector.put(bp, base);
						
					} else {
						base = collector.get(bp);
					}

					Nucleotide nucleotide = Nucleotide.fromCharacter(seq.charAt(j));
					if (nucleotide != null) {
						base.addNucleotidee(nucleotide);
					}
			}
		}
	}

	public Iterator<Base> iterator() {
		return collector.values().iterator();
	}

	/**
	 * Remove bases that don't intercept with the filterRegion.
	 * 
	 * @param filterRegion
	 */
	public void filter(Region filterRegion) {
		
		Iterator<Entry<Long, Base>> iterator = collector.entrySet().iterator();
		
		while (iterator.hasNext()) {
			Entry<Long, Base> entry = iterator.next();
			
			if (!filterRegion.contains(entry.getKey())) {
				iterator.remove();
			}
		}
	}
}
