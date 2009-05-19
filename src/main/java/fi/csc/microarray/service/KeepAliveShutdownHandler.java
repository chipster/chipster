package fi.csc.microarray.service;

import java.util.concurrent.CountDownLatch;

public class KeepAliveShutdownHandler {

	public static void init(final ShutdownCallback shutdownCallback) {
		
		final CountDownLatch keepAliveLatch = new CountDownLatch(1);

		// register shutdown hook
		Runtime.getRuntime().addShutdownHook(new Thread() {
			private CountDownLatch l = keepAliveLatch;

			public void run() {
				l.countDown();
				shutdownCallback.shutdown();
			}
		});
		
		// create and start keep-alive thread
		Thread keepAlive = new Thread() {
			private CountDownLatch l = keepAliveLatch;
			
			public void run() {
				try {
					l.await();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		};
		keepAlive.setDaemon(false);
		keepAlive.start();
	}
}
