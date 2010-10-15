package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.ArrayList;
import java.util.List;

public class Cigar {
	private List<CigarItem> elements = new ArrayList<CigarItem>();
	
	public void addElement(CigarItem e) {
		elements.add(e);
	}
	
	public int getReferenceIndex(int seqIndex) {
		
		int seqCounter = 0;
		int refCounter = 0;
		
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
}
