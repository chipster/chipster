package fi.csc.microarray.messaging;

import java.net.URL;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.UrlListMessage;

/**
 * Reply listener for the url list request.
 * 
 */
public class UrlListMessageListener extends TempTopicMessagingListenerBase {

	private static final Logger logger = Logger.getLogger(UrlListMessageListener.class);

	private List<URL> urlList = null;
	private CountDownLatch latch = new CountDownLatch(1);

	private boolean cancelled = false;

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof UrlListMessage) {
			UrlListMessage urlListMessage = (UrlListMessage) msg;
			this.urlList = urlListMessage.getUrlList();
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return 
	 * @throws AuthCancelledException 
	 * @throws RuntimeException if interrupted
	 */
	public List<URL> waitForReply(long timeout, TimeUnit unit) throws AuthCancelledException {
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
		
		return this.urlList;
	}
	
	@Override
	public void cancel() {
		this.cancelled   = true;
		latch.countDown();
	}
}