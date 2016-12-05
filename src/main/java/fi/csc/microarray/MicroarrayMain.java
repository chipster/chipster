package fi.csc.microarray;

import fi.csc.chipster.ChipsterMain;

/**
 * Wrapper to prevent old scripts etc from breaking.
 * 
 * Use fi.csc.chipster.ChipsterMain instead of this.
 * 
 */
public class MicroarrayMain {

	public static void main(String[] args) {
		System.out.println("This class is obsolete, please use fi.csc.chipster.ChipsterMain instead.");
		ChipsterMain.main(args);
	}
}
