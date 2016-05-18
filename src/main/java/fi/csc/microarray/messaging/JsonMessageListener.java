package fi.csc.microarray.messaging;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.JsonMessage;

public class JsonMessageListener extends TempTopicMessagingListenerBase {

	private JsonMessage value = null;
	private CountDownLatch latch = new CountDownLatch(1);
	private boolean cancelled = false;

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof JsonMessage) {
			JsonMessage jsonMessage = (JsonMessage) msg;
			this.value = jsonMessage;
			latch.countDown();
		}
	}

	/**
	 * @param timeout in given units
	 * @param unit unit of the timeout
	 * @return null if no reply before timeout
	 * @throws AuthCancelledException 
	 * @throws RuntimeException if interrupted
	 */
	public JsonMessage waitForReply(long timeout, TimeUnit unit) throws AuthCancelledException {
		try {
			latch.await(timeout, unit);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} finally {
			// close temp topic
			this.cleanUp();
		}
		if (this.cancelled) {
			throw new AuthCancelledException();
		}
		
		return this.value;
	}

	@Override
	public void cancel() {
		this.cancelled  = true;
		latch.countDown();
	}


}
