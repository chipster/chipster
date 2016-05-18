package fi.csc.microarray.messaging;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fi.csc.microarray.messaging.message.AuthenticationMessage;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public class FirstLoginListener extends TempTopicMessagingListenerBase {

	private CountDownLatch latch = new CountDownLatch(1);
	private boolean cancelled = false;

	public void onChipsterMessage(ChipsterMessage msg) {
		if (msg instanceof AuthenticationMessage && ((AuthenticationMessage)msg).isLoginAck()) {
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
	public void waitForReply() throws AuthCancelledException {
		try {
			latch.await();
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} finally {
			// close temp topic
			this.cleanUp();
		}

		if (this.cancelled) {
			throw new AuthCancelledException();
		}
	}

	@Override
	public void cancel() {
		this.cancelled  = true;
		latch.countDown();
	}

}
