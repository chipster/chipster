package fi.csc.chipster.web.adminweb.ui;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.locks.Lock;

import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.ProgressIndicator;
import com.vaadin.ui.VerticalLayout;

public class AsynchronousView extends VerticalLayout {
	
	private static final int POLLING_INTERVAL = 100;
	
	private Button refreshButton;
	
	private ProgressIndicator progressIndicator = new ProgressIndicator(0.0f);
	private ExecutorService executor = Executors.newCachedThreadPool();
	
	public ProgressIndicator getProggressIndicator() {
		progressIndicator.setWidth(100, Unit.PERCENTAGE);
		return progressIndicator;
	}
	
	protected void waitForUpdate(final Future<?> future, final AfterUpdateCallBack callBack) {				
		
		//This makes the browser start polling, but the browser will get it only if this is executed in this original thread.
		setProgressIndicatorValue(0f);
		
		executor.execute(new Runnable() {
			public void run() {								
				try {									
					/* Separate delay from what happens in the Container, because communication between
					 * threads is messy. Nevertheless, these delays should have approximately same duration
					 * to prevent user from starting several background updates causing concurrent modifications.   
					 */
					final int DELAY = 300; 				
					for (int i = 0; i <= DELAY; i++) {

						try {
							future.get(POLLING_INTERVAL, TimeUnit.MILLISECONDS);
							break;
						} catch (TimeoutException e) {
							//No results yet, update progress bar							
							setProgressIndicatorValue((float)i/DELAY);			
						}
					}
					//Update was successful
					if (callBack != null) {
						callBack.updateDone();
					}

				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				} finally {				
					setProgressIndicatorValue(1.0f);
				}
			}
		});
	}
	
	private void setProgressIndicatorValue(float value) {
		//This happens in initialization 
		if (progressIndicator.getUI() != null ) {
			
			Lock indicatorLock = progressIndicator.getUI().getSession().getLockInstance();
			
			//Component has to be locked before modification from background thread
			indicatorLock.lock();					
			try {
				progressIndicator.setValue(value);
				
				if (value == 1.0f) {
					refreshButton.setEnabled(true);
					progressIndicator.setPollingInterval(Integer.MAX_VALUE);	
				} else {
					refreshButton.setEnabled(false);
					progressIndicator.setPollingInterval(POLLING_INTERVAL);
				}
			} finally {
				indicatorLock.unlock();
			}
		}
	}

	public void submitUpdate(Runnable runnable, AfterUpdateCallBack callBack) {
		Future<?> future = executor.submit(runnable);
		
		waitForUpdate(future, callBack);
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

	public void submitUpdate(Runnable runnable) {
		submitUpdate(runnable, null);
	}
}
