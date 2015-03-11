package fi.csc.microarray.messaging;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

public class SuccessMessageListener extends TempTopicMessagingListenerBase {

	private SuccessMessage message;;
	private CountDownLatch latch = new CountDownLatch(1);

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof SuccessMessage) {
			this.message = (SuccessMessage) msg;
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return null if no reply before timeout
	 * @throws RuntimeException if interrupted
	 */
	public SuccessMessage waitForReply(long timeout, TimeUnit unit) {
		try {
			latch.await(timeout, unit);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} finally {
			// close temp topic
			this.cleanUp();
		}
		return this.message;
	}
}
