package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.analyser.emboss.ACDToSADL;

public class ACDToSADLTest {
    private static String path = "src/test/resources/";

    public static ACDDescription getTestDescription() {
        ACDDescription acd = new ACDDescription();
        acd.fromFile(new File(path + "water.acd"));
        return acd;
    }

    public void run() {
        ACDDescription acd = getTestDescription();
        ACDToSADL converter = new ACDToSADL(acd);
        SADLDescription internal = converter.convert();
        
        // Test number of added parameters
        Assert.assertEquals(internal.parameters().size(), 3);
        Assert.assertEquals(internal.inputs().size(), 3);
        Assert.assertEquals(internal.outputs().size(), 1);
        
        // Test some of the attributes
        Assert.assertEquals(internal.parameters().get(0).getFrom(), "0.0");
        Assert.assertEquals(internal.parameters().get(0).getTo(), "100.0");
    }
    
    public static void main(String[] args) throws Exception {
        new ACDToSADLTest().run();
    }
    
}