package fi.csc.chipster.web.adminweb.ui;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.locks.Lock;

import org.apache.log4j.Logger;

import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.ProgressIndicator;
import com.vaadin.ui.VerticalLayout;

public class AsynchronousView extends VerticalLayout {
	
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(AsynchronousView.class);
	
	private static final int POLLING_INTERVAL = 100;
	
	private Button refreshButton;
	
	private ProgressIndicator progressIndicator = new ProgressIndicator(0.0f);
	private ExecutorService executor = Executors.newCachedThreadPool();

	private long timeout; // seconds
	
	public AsynchronousView(long timeout) {
		this.timeout = timeout;
	}
	
	public AsynchronousView() {
		this(30);
	}

	public ProgressIndicator getProggressIndicator() {
		progressIndicator.setWidth(100, Unit.PERCENTAGE);
		return progressIndicator;
	}
	
	/**
	 * @param future
	 * @param callBack
	 * @param wait Wait always for timeout. Use this when Runnable doesn't block to keep progress indicator refreshing client continuously.
	 */
	protected void waitForUpdate(final Future<?> future, final AfterUpdateCallBack callBack, final boolean wait) {				
		
		//This makes the browser start polling, but the browser will get it only if this is executed in this original thread.
		setProgressIndicatorValue(0f);
		
		executor.execute(new Runnable() {
			public void run() {								
				try {									
					/* Separate delay from what happens in the Container, because communication between
					 * threads is messy. Nevertheless, these delays should have approximately same duration
					 * to prevent user from starting several background updates causing concurrent modifications.   
					 */
					final long POLL_LIMIT = timeout * 1000 / POLLING_INTERVAL;										
					
					for (int i = 0; i <= POLL_LIMIT; i++) {						

						try {
							if (wait) {
								Thread.sleep(POLLING_INTERVAL);
							} else {
								future.get(POLLING_INTERVAL, TimeUnit.MILLISECONDS);
								break;
							}
						} catch (TimeoutException e) {
						}
						//No results yet or wait==true, update progress bar							
						setProgressIndicatorValue((float)i/POLL_LIMIT);
					}
					//Update was successful
					if (callBack != null) {
						updateUI(new Runnable() {
							public void run() {								
								callBack.updateDone();
							}
						});
					}

				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				} finally {				
					setProgressIndicatorValue(1.0f);
				}
			}
		});
	}
	
	private void setProgressIndicatorValue(final float value) {
		
		//Component has to be locked before modification from background thread
		
		this.updateUI(new Runnable() {
			@Override
			public void run() {
				progressIndicator.setValue(value);

				if (value == 1.0f) {
					refreshButton.setEnabled(true);
					progressIndicator.setPollingInterval(Integer.MAX_VALUE);	
				} else {
					refreshButton.setEnabled(false);
					progressIndicator.setPollingInterval(POLLING_INTERVAL);
				}
			}
		});
	}

	
	/**
	 * Thread safe UI updates
	 * 
	 * Execute UI changes in a runnable after attaining the applications's lock.
	 * Runnable won't be run, if the view isn't anymore active. 
	 * 
	 * Make sure this method won't get called before the UI is ready. In 
	 * practice the easiest way to ensure this is to start data queries only 
	 * after the UI is initialized and visible.
	 * 
	 * @param runnable
	 */
	public void updateUI(Runnable runnable) {
		if (this.getUI() == null && this.getParent() == null) {
			// ignore updates to old views, when a new view is already selected
			return;
		}
		
		if (this.getUI().getSession() == null) {
			// ignore updates, when user has already reloaded the page
			return;
		}
		
		Lock lock = this.getUI().getSession().getLockInstance();
		lock.lock();
		
		try {
			runnable.run();
		} finally {
			lock.unlock();
		}
	}

	public void submitUpdate(Runnable runnable, AfterUpdateCallBack callBack, boolean wait) {
		Future<?> future = executor.submit(runnable);
		
		waitForUpdate(future, callBack, wait);
	}
	
	public void submitUpdate(Runnable runnable, AfterUpdateCallBack callBack) {
		submitUpdate(runnable, callBack, false);
	}

	public void submitUpdate(Runnable runnable) {
		submitUpdate(runnable, null, false);
	}
	
	public void submitUpdateAndWait(Runnable runnable) {
		submitUpdate(runnable, null, true);
	}
	
	public Button createRefreshButton(ClickListener listener) {
		refreshButton = new Button("Refresh");
		refreshButton.addClickListener(listener);
		refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
		return refreshButton;			
	}

	public boolean isRefreshButton(Object object) {
		return refreshButton == object;
	}
	
	public long getTimeout() {
		return timeout;
	}
}
