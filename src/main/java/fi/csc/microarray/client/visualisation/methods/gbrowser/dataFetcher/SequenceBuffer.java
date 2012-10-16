package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

public class SequenceBuffer {
	
	private String buffer;
	private long start;
	private Chromosome chr;
	
	public SequenceBuffer (String buffer, long start, Chromosome chr) {
		this.buffer = buffer;
		this.start = start;
		this.chr = chr;
	}
	
	public boolean contains(long requestStart, long requestEnd, Chromosome chr) {
				
		return this.chr.equals(chr) && requestStart >= start && requestEnd <= start + buffer.length();
	}
	
	public String get(long requestStart, long requestEnd, Chromosome chr) {
		if (contains(requestStart, requestEnd, chr)) {
			
			//Create copy of substring to avoid references to original buffer
			String seq = new String(buffer.substring((int)(requestStart - start), (int)(requestEnd - start)));

			return seq;
		}
		return null;
	}
}
