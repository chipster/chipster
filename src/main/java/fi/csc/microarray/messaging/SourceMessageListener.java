package fi.csc.microarray.messaging;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.SourceMessage;

public class SourceMessageListener extends TempTopicMessagingListenerBase {
    
    private CountDownLatch latch = new CountDownLatch(1);
    private String source = null;
    
    
    public String waitForResponse(long timeout, TimeUnit unit) {
        try {
            latch.await(timeout, unit);
            return source;
        } catch (InterruptedException e) {
        	throw new RuntimeException(e);
        }
    }
    
    public void onChipsterMessage(ChipsterMessage message) {
    	
    	// do nothing if already got one reply
    	if (latch.getCount() == 0) {
    		return;
    	}
    	
    	SourceMessage sourceMessage = (SourceMessage) message;
        this.source = sourceMessage.getSource();
        latch.countDown();
    }
}