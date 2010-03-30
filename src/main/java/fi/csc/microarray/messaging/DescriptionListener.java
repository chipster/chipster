package fi.csc.microarray.messaging;

import java.util.Collection;
import java.util.concurrent.CountDownLatch;

import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.DescriptionMessage;

public class DescriptionListener extends TempTopicMessagingListenerBase {
    
    private final CountDownLatch latch = new CountDownLatch(1);
    private Collection<OperationCategory> categories = null;
    private String wantedModule;
    
    public DescriptionListener(String wantedModule) {
        this.wantedModule = wantedModule;
    }
    
    public Collection<OperationCategory> getCategories() {
        return categories;
    }
    
    public void waitForResponse() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    public void onChipsterMessage(ChipsterMessage msg) {
        DescriptionMessage descriptionMsg = (DescriptionMessage) msg;
        
        // TODO change Name to Id
        // TODO check for the right module
        if (descriptionMsg.getModuleName().equals(wantedModule)) {            
            try {
                categories = new OperationGenerator().generateFromMessage(descriptionMsg).values();
            } catch (ParseException e) {
                e.printStackTrace();
            }
            
            latch.countDown();
        }
    }
};