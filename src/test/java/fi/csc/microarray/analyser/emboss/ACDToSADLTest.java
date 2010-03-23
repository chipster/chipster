package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.analyser.emboss.ACDToSADL;
import fi.csc.microarray.analyser.emboss.ACDToSADL.SADLParameterCreator;

public class ACDToSADLTest {
    private static String path = "src/test/resources/";

    public static ACDDescription getTestDescription(String appName) {
        ACDDescription acd = new ACDDescription(new File(path + appName + ".acd"));
        return acd;
    }

    @Test
    public void testACDToSADL() {
        ACDDescription acd = getTestDescription("water");
        ACDToSADL converter = new ACDToSADL(acd);
        SADLDescription internal = converter.convert();
        
        // Test number of added parameters
        Assert.assertEquals(internal.parameters().size(), 3);
        Assert.assertEquals(internal.inputs().size(), 3);
        Assert.assertEquals(internal.outputs().size(), 1);
        
        // Test some of the attributes
        Assert.assertEquals(internal.parameters().get(0).getFrom(), "0.0");
        Assert.assertEquals(internal.parameters().get(0).getTo(), "100.0");
        
        // Test parameter creation
        ACDParameter acdParam = new ACDParameter("list", "param", "", null);
        String[] titles = {"all", "your", "base"};
        String[] values = {"a", "y", "b"};
        acdParam.setList(titles, values);
        Parameter sadlParam = SADLParameterCreator.createParameter(acdParam);
        Assert.assertEquals(sadlParam.getSelectionOptions()[0], "all"); 
    }
    
    public static void main(String[] args) throws Exception {
        new ACDToSADLTest().testACDToSADL();
    }
    
}