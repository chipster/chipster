package fi.csc.microarray.util;

import java.io.IOException;

import de.schlichtherle.truezip.zip.ZipFile;

/**
 * Separate from IOUtils to avoid truezip dependency there.
 * 
 * @author klemela
 */
public class ZipUtils {
	public static void closeIfPossible(ZipFile zipFile) {
		if (zipFile != null) {
			try {
				zipFile.close();
			} catch (IOException e) {
				// ignore
			}
		}
	}
}
