package fi.csc.microarray.util;

import net.sf.samtools.SAMRecord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserSettings.CoverageType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

public class BamUtils {

	public static Strand getStrand(SAMRecord record, CoverageType coverageType) {
		if (coverageType == CoverageType.STRAND) {
			
			if (record.getReadNegativeStrandFlag()) {
				return Strand.REVERSE;
			} else {
				return Strand.FORWARD;
			}
			
		} else if (coverageType == CoverageType.STRAND_XS) {
			
			char xs = record.getCharacterAttribute("XS");
			if (xs == '+') {
				return Strand.FORWARD;
			} else if (xs == '-') {
				return Strand.REVERSE;
			}
		}
		return null;
	}
}
