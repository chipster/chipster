package fi.csc.microarray.filebroker;

import java.io.IOException;
import java.io.InputStream;
import java.net.URLConnection;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Timer;
import java.util.TimerTask;

import javax.xml.bind.DatatypeConverter;

/**
 * <p>Calculate md5 checksum of the stream to verify its contents. This class provides a transparent
 * input stream that calculates a checksum of the data passing through the input stream.</p>
 * 
 * <p>It's not practical to read through a big file just to calculate checksum, but it's quite cheap to calculate 
 * md5 on the fly when the file is read anyway for some other purpose. Still the calculation consumes some resources and 
 * therefore the constructors have a parameter useChecksums that can be used to bypass checksum calculation 
 * to attain higher throughput (checksum calculation will top out around 200 MB/s utilizing fully one CPU core).</p>
 * 
 * This class can be used in three different use cases:
 * <ul>
 * <li>Local calculation</li>
 * <li>Download</li>
 * <li>Upload</li>
 * </ul>
 * 
 * <h2>Local calculation</h2>
 * 
 * Use two parameter constructor to create an instance of this class, read through the stream and call getChecksum()
 * to get the md5 as String. See method Md5FileUtils.calculateMd5() if you don't have any other use for the stream, but 
 * just want to calculate the checksum.
 * 
 * <h2>Download</h2>
 * 
 * This class can do the comparison of checksums if you are downloading a file over HTTP 
 * and the server sends the original checksum in 'Etag' HTTP header. Use three parameter constructor
 * that takes an URLConnection object and then read through the stream. By calling verifyChecksums(), the 
 * locally calculated checksum is compared to checksum sent by the server.  
 * 
 * <h2>Upload</h2>
 * 
 * There is no separate OutputStream implementation of this class, but if you are reading the data 
 * from InputStream, you can still calculate the checksum from that InputStream and compare it with
 * the server response of your outbound HTTP connection.
 * 
 * @author klemela
 *
 */
public class ChecksumInputStream extends DigestInputStream {
	
	public static final String HTTP_CHECKSUM_KEY = "Etag";
	
	private boolean useChecksums;
	private URLConnection connection;
	private String checksum = null;
	private long bytes = 0;

	private boolean isClosed;
	
	public ChecksumInputStream(InputStream baseStream, boolean useChecksums, final URLConnection connection) {
		super(baseStream, getDigestWithRuntimeException());
		
		this.useChecksums = useChecksums;
		super.on(useChecksums);
		
		this.connection = connection;
		
		/*
		 * Debug unclosed streams
		 * 
		 * We often read only the beginning of the file, so it's important to
		 * close the stream to avoid leaving threads waiting on the server.
		 * checkClose() waits for 5 seconds and prints a stack trace if the
		 * stream isn't closed. Obviously this will report false positives if
		 * the files are big. More sophisticated check would track the stream
		 * usage and issue a warning when the stream is left idle. This would
		 * allow enabling this in production.
		 */
		//checkClose(Exceptions.getStackTrace(Thread.currentThread().getStackTrace()));
	}

	@SuppressWarnings("unused")
	private void checkClose(final String stack) {
		Timer timer = new Timer();
		timer.schedule(new TimerTask() {
			@Override
			public void run() {
				if (!isClosed) {
					System.out.println("Unclosed input stream found");
					System.out.println(stack);
					System.out.println();
				}
			}
		}, 5_000);
	}

	/**
	 * There is not much clients can do in case of NoSuchAlgorithmException, so there is no point to have it
	 * in method signatures. Convert it to RuntimeException here because its not possible to make try-catch block
	 * around the super constructor call.
	 * 
	 * @return
	 */
	private static MessageDigest getDigestWithRuntimeException() {
		try {
			return MessageDigest.getInstance("MD5");
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
	}

	public ChecksumInputStream(InputStream baseStream, boolean useChecksums) {
		this(baseStream, useChecksums, null);
	}

	@Override
	public int read() throws IOException {		
		checkState();		
		int b = super.read();
		if (b != -1) { //if not end of stream
			bytes += 1; //one byte read
		}
		return b;		
	}
	
	@Override
	public int read(byte b[], int off, int len) throws IOException {
		checkState();
		int bytesRead = super.read(b, off, len);
		if (bytesRead != -1) {
			bytes += bytesRead;
		}
		return bytesRead;
	}
		
	@Override
	public int available() throws IOException {
		checkState();
		return super.available();
	}
		
	private void checkState() {
		if (useChecksums && checksum != null) {
			throw new IllegalStateException("use of InputStream is prevented after checksum calculation");
		}
	}

	/**
	 * This can be called only once because this resets the digest. Subsequent calls will throw IllegalStateException.
	 * 
	 * @return
	 */
	public String getChecksum() {
		checkState();
		
		if (useChecksums) {
			checksum =  DatatypeConverter.printHexBinary(super.getMessageDigest().digest()).toLowerCase();
			return checksum;
		} else {
			return null;
		}
	}
		
	private String getRemoteChecksum() {
		if (connection == null) {
			throw new IllegalStateException("URLConnection wasn't given when this object was constructed: server checksum is unavailable");
		}
		return connection.getHeaderField(HTTP_CHECKSUM_KEY);		
	}
	
	/**
	 * Verify that checksums in both ends of network transmission are equal.
	 * 
	 * @return locally calculated checksum or null id calculation is disabled
	 * @throws ChecksumException
	 */
	public String verifyChecksums() throws ChecksumException {
		
		if (useChecksums) {
			String localChecksum = getChecksum();
			Md5FileUtils.verify(getRemoteChecksum(), localChecksum);
			return localChecksum;
		}
		return null;
	}

	/**
	 * 
	 * Verifies that contentLength equals to the stream length, when non-null value is given. This is a
	 * backup verification method, because checksums can be disabled.
	 * 
	 * @throws ContentLengthException 
	 */
	public void verifyContentLength(Long contentLength) throws ChecksumException, ContentLengthException {				
		
		Md5FileUtils.verify(bytes, contentLength);
	}

	/**
	 * For normal use verifyContentLength() should be enough, but this can 
	 * be useful for debug messages.
	 * 
	 * @return
	 */
	public Long getContentLength() {
		return bytes;
	}
	
	@Override
	public void close() throws IOException {
		this.isClosed = true;
		super.close();
	}
}
