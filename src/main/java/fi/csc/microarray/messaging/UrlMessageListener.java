package fi.csc.microarray.messaging;

import java.net.URL;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.UrlMessage;

/**
 * Reply listener for the url request.
 * 
 */
public class UrlMessageListener extends TempTopicMessagingListenerBase {

	private static final Logger logger = Logger.getLogger(UrlMessageListener.class);

	private URL newUrl = null;
	private CountDownLatch latch = new CountDownLatch(1);

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof UrlMessage) {
			UrlMessage urlMessage = (UrlMessage) msg;
			this.newUrl = urlMessage.getUrl();
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return may be null
	 * @throws RuntimeException if interrupted
	 */
	public URL waitForReply(long timeout, TimeUnit unit) {
		try {
			latch.await(timeout, unit);
		} catch (InterruptedException e) {
			logger.warn("interrupted while waiting for latch", e);
		}
		return this.newUrl;
	}
}