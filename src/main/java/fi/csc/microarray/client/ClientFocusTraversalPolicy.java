package fi.csc.microarray.client;

import java.awt.Component;
import java.awt.Container;
import java.awt.FocusTraversalPolicy;
import java.util.Vector;

import org.apache.log4j.Logger;

public class ClientFocusTraversalPolicy extends FocusTraversalPolicy{
    Vector<Component> order;
    
    private static final Logger logger = Logger.getLogger(ClientFocusTraversalPolicy.class); 

    public ClientFocusTraversalPolicy(Vector<Component> order) {
        this.order = new Vector<Component>(order.size());
        this.order.addAll(order);
    }
    public Component getComponentAfter(Container focusCycleRoot,
                                       Component aComponent)
    {
        int idx = (order.indexOf(aComponent) + 1) % order.size();
        logger.debug("Focus moved to next component: " + order.get(idx));
        Component comp = order.get(idx);
        if(!comp.isEnabled()){
        	comp = getComponentAfter(focusCycleRoot, comp);
        }
        return comp;        
    }

    public Component getComponentBefore(Container focusCycleRoot,
                                        Component aComponent)
    {
        int idx = order.indexOf(aComponent) - 1;
        if (idx < 0) {
            idx = order.size() - 1;
        }
        logger.debug("Focus moved to previous component: " + order.get(idx));
        Component comp = order.get(idx);
        if(!comp.isEnabled()){
        	comp = getComponentBefore(focusCycleRoot, comp);
        }
        return comp;
    }

    public Component getDefaultComponent(Container focusCycleRoot) {
    	logger.debug("Focus moved to default component: " + order.get(0));
    	Component comp = order.get(0);
        if(!comp.isEnabled()){
        	comp = getComponentAfter(focusCycleRoot, comp);
        }
        return comp;
    }

    public Component getLastComponent(Container focusCycleRoot) {
    	logger.debug("Focus moved to last component: " + order.lastElement());
        Component comp = order.lastElement();
        if(!comp.isEnabled()){
        	comp = getComponentBefore(focusCycleRoot, comp);
        }
        return comp;
    }

    public Component getFirstComponent(Container focusCycleRoot) {
    	logger.debug("Focus moved to first component: " + order.get(0));
        return getDefaultComponent(focusCycleRoot);
    }
}
