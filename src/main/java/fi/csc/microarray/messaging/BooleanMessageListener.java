package fi.csc.microarray.messaging;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fi.csc.microarray.messaging.message.BooleanMessage;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public class BooleanMessageListener extends TempTopicMessagingListenerBase {

	private Boolean value = null;
	private CountDownLatch latch = new CountDownLatch(1);

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof BooleanMessage) {
			BooleanMessage booleanMessage = (BooleanMessage) msg;
			this.value = booleanMessage.getValue();
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return null if no reply before timeout
	 * @throws RuntimeException if interrupted
	 */
	public Boolean waitForReply(long timeout, TimeUnit unit) {
		try {
			latch.await(timeout, unit);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} finally {
			// close temp topic
			this.cleanUp();
		}
		return this.value;
	}
}
