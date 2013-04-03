package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

/**
 * DNA strand, either forward or reverse. Value BOTH is usually used in tracks to denote that
 * information from both strands should be included in a single track.
 * 
 * @author Petri Klemelä
 */
public enum Strand {
	FORWARD, 
	REVERSE, 
	BOTH,
	NONE,
	UNRECOGNIZED;
	
	@Override
	public String toString() {
		switch (this) {
		case FORWARD:
			return "+";
		case REVERSE:
			return "-";
		default:
			return super.toString();
		}
	}
}
