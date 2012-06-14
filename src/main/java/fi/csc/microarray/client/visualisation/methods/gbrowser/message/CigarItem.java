package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;

/**
 * Single item (continuous block of same type bases) of CIGAR string.
 * 
 * @author Aleksi Kallio
 *
 */
public class CigarItem {

	public enum CigarItemType {
		D, //Deletion
		EQ, //Matches the reference
		H, //Hard clip
		I, //Insertion vs.
		M, //Match or mismatch
		N, //Skipped region from the reference.
		P, //Padding.
		S, //Soft clip.
		X } //Mismatches the reference.

	private CigarElement cigarElement;

	public CigarItem(CigarElement cigarElement) {
		this.cigarElement = cigarElement;
	}

	public long getLength() {
		return cigarElement.getLength();
	}

	public String getType() {
		return cigarElement.getOperator().toString();
	}

	public boolean isVisible() {
		return cigarElement.getOperator() == CigarOperator.M || cigarElement.getOperator() == CigarOperator.X || cigarElement.getOperator() == CigarOperator.EQ;
	}

	public boolean consumesReferenceBases() {
		return cigarElement.getOperator().consumesReferenceBases();
	}

	public boolean consumesReadBases() {
		return cigarElement.getOperator().consumesReadBases();
	}

	public CigarItemType getCigarItemType() {
		//Translate sam-library specific type into generic Chipster type
		switch (cigarElement.getOperator()) {
		case D:
			return CigarItemType.D;
		case EQ:
			return CigarItemType.EQ;
		case H:
			return CigarItemType.H;
		case I:
			return CigarItemType.I;
		case M:
			return CigarItemType.M;
		case N:
			return CigarItemType.N;
		case P:
			return CigarItemType.P;
		case S:
			return CigarItemType.S;
		case X:
			return CigarItemType.X;

		default:
			return null;
		}
	}
}
