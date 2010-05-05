package fi.csc.microarray.messaging;

import java.util.Collection;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.DescriptionMessage;

public class DescriptionMessageListener extends TempTopicMessagingListenerBase {
    
	private static final Logger logger = Logger.getLogger(DescriptionMessageListener.class);
    
	private final CountDownLatch latch = new CountDownLatch(1);
    private Collection<OperationCategory> categories = null;
    private String wantedModule;
    private boolean finished = false;
    
    public DescriptionMessageListener(String wantedModule) {
        this.wantedModule = wantedModule;
    }
    
    public Collection<OperationCategory> getCategories() {
        return categories;
    }
    
    public void waitForResponse() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            logger.warn("interrupted while waiting for descriptions message", e);
        }
    }
    
    public void onChipsterMessage(ChipsterMessage msg) {
        if (finished) {
        	return;
        }
    	
    	DescriptionMessage descriptionMsg = (DescriptionMessage) msg;
        
        // TODO change Name to Id
        // TODO check for the right module
    	// only wait for one description message for the correct module
    	if (descriptionMsg.getModuleName().equals(wantedModule)) {            
        	try {
                categories = new OperationGenerator().generateFromMessage(descriptionMsg).values();
            } catch (ParseException e) {
                logger.warn("parsing descriptions message failed", e);
            } finally {
                finished = true;
                latch.countDown();
            }
        }
    }
};