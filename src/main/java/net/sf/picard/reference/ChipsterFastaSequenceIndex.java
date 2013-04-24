/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package net.sf.picard.reference;

import java.io.IOException;

import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.FastaSequenceIndexEntry;
import net.sf.samtools.SAMSequenceRecord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;

/**
 * Copy of Picard class modified to support urls 
 */
public class ChipsterFastaSequenceIndex extends FastaSequenceIndex implements Iterable<FastaSequenceIndexEntry> {

	public ChipsterFastaSequenceIndex( LineDataSource indexFile ) {
        parseIndexFile(indexFile);
    }

	private void parseIndexFile(LineDataSource indexFile) {
        try {

            int sequenceIndex = 0;

            String line = null;
            while ((line = indexFile.readLine()) != null) {

            	String[] splitted = line.split("\t");
                // Parse the index line.
                String contig = splitted[0];
                long size = Long.valueOf(splitted[1]);
                long location = Long.valueOf(splitted[2]);
                int basesPerLine = Integer.valueOf(splitted[3]);
                int bytesPerLine = Integer.valueOf(splitted[4]);

                contig = SAMSequenceRecord.truncateSequenceName(contig);
                // Build sequence structure
                add(new FastaSequenceIndexEntry(contig,location,size,basesPerLine,bytesPerLine, sequenceIndex++) );
            }
        } catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
}
