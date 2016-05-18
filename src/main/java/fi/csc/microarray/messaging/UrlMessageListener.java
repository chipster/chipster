package fi.csc.microarray.messaging;

import java.net.URL;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.filebroker.QuotaExceededException;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;
import fi.csc.microarray.messaging.message.UrlMessage;

/**
 * Reply listener for the url request.
 * 
 */
public class UrlMessageListener extends TempTopicMessagingListenerBase {

	private static final Logger logger = Logger.getLogger(UrlMessageListener.class);

	private URL newUrl = null;
	private FileBrokerException exception;
	private CountDownLatch latch = new CountDownLatch(1);

	private boolean cancelled = false;

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof UrlMessage) {
			UrlMessage urlMessage = (UrlMessage) msg;
			this.newUrl = urlMessage.getUrl();
			latch.countDown();
		}
		if (msg instanceof SuccessMessage) {
			SuccessMessage successMessage = (SuccessMessage) msg;
			if (FileServer.ERROR_QUOTA_EXCEEDED.equals(successMessage.getErrorMessage())) {
				exception = new QuotaExceededException();
			}
			latch.countDown();
		}
		if (msg instanceof CommandMessage) {
			CommandMessage commandMessage = (CommandMessage) msg;
			// the following is most likely CommandMessage.COMMAND_FILE_OPERATION_DENIED
			exception = new FileBrokerException(commandMessage.getCommand());  
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return may be null
	 * @throws FileBrokerException if file operation failed
	 * @throws AuthCancelledException 
	 * @throws RuntimeException if interrupted
	 */
	public URL waitForReply(long timeout, TimeUnit unit) throws FileBrokerException, AuthCancelledException {
		try {
			latch.await(timeout, unit);
		} catch (InterruptedException e) {
			logger.warn("interrupted while waiting for latch", e);
		} finally {
			// close temp topic
			this.cleanUp();
		}
		
		if (this.cancelled) {
			throw new AuthCancelledException();
		}
		
		if (exception != null) {
			throw exception;
		}
		
		return this.newUrl;
	}
	
	@Override
	public void cancel() {
		this.cancelled   = true;
		latch.countDown();
	}
}