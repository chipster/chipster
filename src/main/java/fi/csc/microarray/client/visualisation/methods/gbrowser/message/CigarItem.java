package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;

public class CigarItem {
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
		return cigarElement.getOperator() == CigarOperator.M;
	}
	
	public boolean consumesReferenceBases() {
		return cigarElement.getOperator().consumesReferenceBases();
	}

	public boolean consumesReadBases() {
		return cigarElement.getOperator().consumesReadBases();
	}
}
