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

package htsjdk.samtools.reference;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ByteDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.LineDataSource;
import htsjdk.samtools.SAMException;

/**
 * Copy of Picard class modified to support urls and 
 *  - smaller buffer size to avoid downloading useless extra bytes)
 *  - access to FastaSequenceIndex
 */
public class ChipsterIndexedFastaSequenceFile extends PicardIndexedFastaSequenceFile {
 
	private final ByteDataSource channel;

    private final htsjdk.samtools.reference.FastaSequenceIndex index;
    
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

    /**
     * Gets the subsequence of the contig in the range [start,stop]
     * @param contig Contig whose subsequence to retrieve.
     * @param start inclusive, 1-based start of region.
     * @param stop inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
        if(start > stop + 1)
            throw new SAMException(String.format("Malformed query; start point %d lies after end point %d",start,stop));

        FastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

        if(stop > indexEntry.getSize())
            throw new SAMException("Query asks for data past end of contig");

        int length = (int)(stop - start + 1);

        byte[] target = new byte[length];
        ByteBuffer targetBuffer = ByteBuffer.wrap(target);

        final int basesPerLine = indexEntry.getBasesPerLine();
        final int bytesPerLine = indexEntry.getBytesPerLine();
        final int terminatorLength = bytesPerLine - basesPerLine;

        long startOffset = ((start-1)/basesPerLine)*bytesPerLine + (start-1)%basesPerLine;
        // Cast to long so the second argument cannot overflow a signed integer.
        final long minBufferSize = Math.min((long) BUFFER_SIZE, (long)(length / basesPerLine + 2) * (long)bytesPerLine);
        if (minBufferSize > Integer.MAX_VALUE) throw new SAMException("Buffer is too large: " +  minBufferSize);

        // Allocate a buffer for reading in sequence data.
        final ByteBuffer channelBuffer = ByteBuffer.allocate((int)minBufferSize);

        while(targetBuffer.position() < length) {
            // If the bufferOffset is currently within the eol characters in the string, push the bufferOffset forward to the next printable character.
            startOffset += Math.max((int)(startOffset%bytesPerLine - basesPerLine + 1),0);

            try {
                //startOffset += readFromPosition(channel, channelBuffer, indexEntry.getLocation()+startOffset);
            	byte[] bytes = channel.read(indexEntry.getLocation() + startOffset, channelBuffer.remaining());
            	channelBuffer.put(bytes);
            	startOffset += bytes.length;
            }
            catch(IOException ex) {
                throw new SAMException("Unable to load " + contig + "(" + start + ", " + stop + ") from " + getAbsolutePath(), ex);
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
