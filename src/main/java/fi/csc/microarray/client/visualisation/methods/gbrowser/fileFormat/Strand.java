package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

/**
 * DNA strand.
 * 
 * Value BOTH is usually used in tracks to denote that
 * information from both strands should be included in
 * a single track.
 */
public enum Strand {
	FORWARD, 
	REVERSED, 
	BOTH;
	
	@Override
	public String toString() {
		switch (this) {
		case FORWARD:
			return "+";
		case REVERSED:
			return "-";
		default:
			return super.toString();
		}
	}
}
