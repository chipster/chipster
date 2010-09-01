package fi.csc.microarray.util;

import org.apache.log4j.Logger;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

public class StreamStartCache implements InputStreamSource {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(StreamStartCache.class);
	
	private static class StartCachedInputStream extends InputStream {

		private ByteArrayInputStream cachedStream;
		private InputStream stream = null;
		private InputStreamSource source;	
		private boolean useCache = true;
		private int cachePosition = 0;
		
		StartCachedInputStream(ByteArrayInputStream cachedStream, InputStreamSource source) {
			this.cachedStream = cachedStream;
			this.source = source;
		}
		
		@Override
		public int read() throws IOException {
			
			// check if cache has become unusable 
			if (useCache && cachedStream.available() == 0) {
				
				logger.debug("stream start cache miss");
				useCache = false;
				
				// initialise actual stream
				this.stream = source.getInputStream();
				
				// go to correct spot
				for (int i = 0; i < cachePosition; i++) {
					stream.read();
				}
			}
			
			
			if (useCache) {
				cachePosition++;
				return cachedStream.read();
			} else {
				return stream.read();
			}
		}
	}
	
	public static final int CACHE_SIZE = 10 * 1024; // 10k

	private byte[] cache;
	private InputStreamSource source;
	
	public StreamStartCache(InputStream initialiser, InputStreamSource source) throws IOException {
		byte[] tempCache = new byte[CACHE_SIZE];
		int read = initialiser.read(tempCache);
		if (read == -1) {
			read = 0; // if EOF, zero bytes were read
		}
		this.cache = new byte[read];
		System.arraycopy(tempCache, 0, cache, 0, read);
		
		this.source = source;
	}
	
	private ByteArrayInputStream getByteArrayInputStream() {
		return new ByteArrayInputStream(cache);
	}

	public InputStream getInputStream() {
		return new StartCachedInputStream(getByteArrayInputStream(), source);
	}
}
