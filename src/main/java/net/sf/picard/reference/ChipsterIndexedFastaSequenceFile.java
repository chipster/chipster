/*
 * 
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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.PicardException;
import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.FastaSequenceIndexEntry;
import net.sf.picard.reference.ReferenceSequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ByteDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;

/**
 * Copy of Picard class modified to support urls and 
 *  - smaller buffer size to avoid downloading useless extra bytes)
 *  - access to FastaSequenceIndex
 */
public class ChipsterIndexedFastaSequenceFile extends PicardIndexedFastaSequenceFile {
 
	private final ByteDataSource channel;

    private final FastaSequenceIndex index;
    
    /**
     * Size of the read buffer.
     */
    protected static final int BUFFER_SIZE = 1 * 1024;

     public ChipsterIndexedFastaSequenceFile(final ByteDataSource file, final FastaSequenceIndex index) {

        super();
        this.index = index;
        channel = file;
    }

    public ChipsterIndexedFastaSequenceFile(ByteDataSource data, LineDataSource index) throws FileNotFoundException {
        this(data, new ChipsterFastaSequenceIndex(index));
    }

    public List<String> getContigs() {
    	
    	List<String> contigs = new LinkedList<String>();
    	
    	for (FastaSequenceIndexEntry entry : index) {
    		contigs.add(entry.getContig());
    	}
    	return contigs;
    }

    public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
        if(start > stop + 1)
            throw new PicardException(String.format("Malformed query; start point %d lies after end point %d",start,stop));

        FastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

        if(stop > indexEntry.getSize())
            throw new PicardException("Query asks for data past end of contig");

        int length = (int)(stop - start + 1);

        byte[] target = new byte[length];
        ByteBuffer targetBuffer = ByteBuffer.wrap(target);

        final int basesPerLine = indexEntry.getBasesPerLine();
        final int bytesPerLine = indexEntry.getBytesPerLine();
        final int terminatorLength = bytesPerLine - basesPerLine;

        long startOffset = ((start-1)/basesPerLine)*bytesPerLine + (start-1)%basesPerLine;

        // Allocate a 128K buffer for reading in sequence data.
        ByteBuffer channelBuffer = ByteBuffer.allocate(BUFFER_SIZE);

        while(targetBuffer.position() < length) {
            // If the bufferOffset is currently within the eol characters in the string, push the bufferOffset forward to the next printable character.
            startOffset += Math.max((int)(startOffset%bytesPerLine - basesPerLine + 1),0);

            try {
                 //startOffset += channel.read(channelBuffer,indexEntry.getLocation()+startOffset);
            	startOffset += channel.read(indexEntry.getLocation()+startOffset, channelBuffer.array());
            }
            catch(IOException ex) {
                throw new PicardException("Unable to load " + contig + "(" + start + ", " + stop + ") from " + channel.toString(), ex);
            }

            // Reset the buffer for outbound transfers.
            channelBuffer.flip();

            // Calculate the size of the next run of bases based on the contents we've already retrieved.
            final int positionInContig = (int)start-1+targetBuffer.position();
            final int nextBaseSpan = Math.min(basesPerLine-positionInContig%basesPerLine,length-targetBuffer.position());
            // Cap the bytes to transfer by limiting the nextBaseSpan to the size of the channel buffer.
            int bytesToTransfer = Math.min(nextBaseSpan,channelBuffer.capacity());

            channelBuffer.limit(channelBuffer.position()+bytesToTransfer);

            while(channelBuffer.hasRemaining()) {
                targetBuffer.put(channelBuffer);

                bytesToTransfer = Math.min(basesPerLine,length-targetBuffer.position());
                channelBuffer.limit(Math.min(channelBuffer.position()+bytesToTransfer+terminatorLength,channelBuffer.capacity()));
                channelBuffer.position(Math.min(channelBuffer.position()+terminatorLength,channelBuffer.capacity()));
            }

            // Reset the buffer for inbound transfers.
            channelBuffer.flip();
        }

        return new ReferenceSequence( contig, indexEntry.getSequenceIndex(), target );
    }
}
