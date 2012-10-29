package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.emboss.ACDToSADL.SADLParameterCreator;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLParser.ParseException;

public class ACDToSADLTest {
    private static String path = "src/test/resources/";
    
    @BeforeTest
    protected void setUp() throws Exception {
    	DirectoryLayout.uninitialise();
        DirectoryLayout.initialiseSimpleLayout();
    }

    public static ACDDescription getTestDescription(String appName) {
        ACDDescription acd = new ACDDescription(new File(path + appName + ".acd"));
        return acd;
    }

    @Test
    public void testACDToSADL() {
        // Load water.acd
        ACDDescription acd = getTestDescription("water");
        SADLDescription sadl = ACDToSADL.convert(acd, "water.acd");
        
        // Test number of added parameters
        Assert.assertEquals(sadl.parameters().size(), 3);
        Assert.assertEquals(sadl.inputs().size(), 2);
        Assert.assertEquals(sadl.outputs().size(), 1);
        
        // Test some of the attributes
        Assert.assertEquals(sadl.parameters().get(0).getFrom(), "0.0");
        Assert.assertEquals(sadl.parameters().get(0).getTo(), "100.0");
        
        // Test parameter creation
        ACDParameter acdParam = new ACDParameter("list", "param", "", null);
        String[] titles = {"all", "your", "base"};
        String[] values = {"a", "y", "b"};
        acdParam.setList(titles, values);
        acdParam.setAttribute("information", "List param");
        acdParam.setAttribute("help", "Long description");
        Parameter sadlParam = SADLParameterCreator.createParameter(acdParam);
        Assert.assertEquals(sadlParam.getName().getDisplayName(), "List param");
        Assert.assertEquals(sadlParam.getComment(), "Long description");
        Assert.assertEquals(sadlParam.getSelectionOptions()[0].toString(),
                            Name.createName("a", "all").toString());
        
        // Load complex.acd
        acd = getTestDescription("complex");
        sadl = ACDToSADL.convert(acd, "complex.acd");
        
        // Test multiple defaults in selection lists
        Assert.assertEquals(acd.getParameter("multiple").getDefaults().length, 3);
        for (Parameter parameter : sadl.parameters()) {
            if (parameter.getName().getID().equals("multiple")) {
                Assert.assertEquals(parameter.getFrom(), "1");
                Assert.assertEquals(parameter.getTo(), "4");
                Assert.assertEquals(parameter.getDefaultValues().length, 3);
                Assert.assertEquals(parameter.getDefaultValues()[0], "b");
                Assert.assertEquals(parameter.getDefaultValues()[1], "o");
                Assert.assertEquals(parameter.getDefaultValues()[2], "g");
            }
        }
        
        // Try converting to string and parsing again
        try {
            sadl = new SADLParser().parse(sadl.toString());
            for (Parameter parameter : sadl.parameters()) {
                if (parameter.getName().getID().equals("multiple")) {
                    Assert.assertEquals(parameter.getFrom(), "1");
                    Assert.assertEquals(parameter.getTo(), "4");
                    Assert.assertEquals(parameter.getDefaultValues().length, 3);
                    Assert.assertEquals(parameter.getDefaultValues()[0], "b");
                    Assert.assertEquals(parameter.getDefaultValues()[1], "o");
                    Assert.assertEquals(parameter.getDefaultValues()[2], "g");
                }
            }
        } catch (ParseException e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) throws Exception {
        ACDToSADLTest test = new ACDToSADLTest();
        test.setUp();
        test.testACDToSADL();
    }
    
}