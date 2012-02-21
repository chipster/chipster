package fi.csc.microarray.analyser.emboss;

import java.io.File;
import java.util.LinkedHashMap;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.emboss.ACDDescription;

public class ACDDescriptionTest {
    ACDDescription acd;
    
    String path = "src/test/resources/";
    
    @BeforeClass
    public void setUp() {
        acd = new ACDDescription(new File(path + "emma.acd"));
    }
    
    @Test
    public void testACDDescription() {
        
        // Test application description
        Assert.assertEquals(acd.getName(), "emma");
        Assert.assertEquals(acd.getDescription(), "Multiple sequence alignment (ClustalW wrapper)");
        
        // Test application groups
        Assert.assertEquals(acd.getGroups().size(), 1);
        Assert.assertEquals(acd.getGroups().get(0), "Alignment:Multiple");
        
        // Test parameter sections
        Assert.assertEquals(acd.getParameters("input", null, false).size(), 5);
        Assert.assertEquals(acd.getParameters("input", null, true).size(), 11);
        Assert.assertEquals(acd.getParameters("input", "matrixsection", false).size(), 3);
        Assert.assertEquals(acd.getParameters("additional", null, true).size(), 15);
        
        // Test parameter attributes
        Assert.assertEquals(acd.getParameter("maxdiv").getAttribute("default"), "30");
        Assert.assertEquals(acd.getParameter("dendoutfile").getAttribute("extension"), "dnd");
        Assert.assertEquals(acd.getParameter("dnamatrix").getList().size(), 3);
        
        // Test parameter obligatoriness
        Assert.assertFalse(acd.getParameter("maxdiv").isRequired());
        Assert.assertTrue(acd.getParameter("maxdiv").isAdditional());
        Assert.assertTrue(acd.getParameter("sequence").isRequired());
        Assert.assertFalse(acd.getParameter("sequence").isAdditional());
        Assert.assertFalse(acd.getParameter("sequence").isAdvanced());
        
        // Test parameter validation using acd file
        Assert.assertTrue(acd.getParameter("pwgapopen").validate("0.1"));
        Assert.assertFalse(acd.getParameter("pwgapopen").validate("-0.1"));
    }
    
    @Test
    public void testACDParameter() {
        ACDParameter param;
        
        // Expression evaluation
        param = new ACDParameter("integer", "param", "", null);
        LinkedHashMap<String, String> map = new LinkedHashMap<String, String>();
        map.put("foo", "n");
        map.put("bar", "21");
        Assert.assertEquals(ACDParameter.resolveExp("$(foo)", map), "false");
        Assert.assertEquals(ACDParameter.resolveExp("@(!$(foo))", map), "true");
        Assert.assertEquals(ACDParameter.resolveExp("@($(bar)+2)", map), "23");
        
        // Value normalization
        param = new ACDParameter("boolean", "param", "", null);
        Assert.assertEquals(param.normalize("yes"), "Y");
        Assert.assertEquals(param.normalize("n"), "N");
        Assert.assertFalse(param.normalize("taip").equals("Y"));
        param = new ACDParameter("selection", "param", "", null);
        Assert.assertEquals(param.normalize("1,2,3"), "1;2;3");
        param.setAttribute("delimiter", "-");
        Assert.assertEquals(param.normalize("1,2,3"), "1-2-3");

        // Array validation
        param = new ACDParameter("array", "param", "", null);
        Assert.assertTrue(param.validate("1.5,1.0,0.5"));
        Assert.assertTrue(param.validate("1.5 1.0 2"));
        Assert.assertFalse(param.validate("1.5 1.0 a"));
        
        // Boolean validation
        param = new ACDParameter("boolean", "param", "", null);
        Assert.assertTrue(param.validate("true"));
        Assert.assertTrue(param.validate("Y"));
        Assert.assertTrue(param.validate("0"));
        Assert.assertFalse(param.validate("X"));
        Assert.assertFalse(param.validate("24"));
        
        // Float validation
        param = new ACDParameter("float", "param", "", null);
        param.setAttribute("minimum", "-2.2");
        param.setAttribute("maximum", "2");
        Assert.assertTrue(param.validate("-2.2"));
        Assert.assertFalse(param.validate("-2.3"));
        Assert.assertFalse(param.validate("2.01"));
        
        // Integer validation
        param = new ACDParameter("integer", "param", "", null);
        param.setAttribute("minimum", "-7");
        param.setAttribute("maximum", "23");
        Assert.assertTrue(param.validate("23"));
        Assert.assertTrue(param.validate("-7"));
        Assert.assertFalse(param.validate("-2.3"));
        Assert.assertFalse(param.validate("24"));
        
        // List validation
        param = acd.getParameter("pwmatrix");
        Assert.assertTrue(param.validate("o"));
        Assert.assertFalse(param.validate("x"));
        param = new ACDParameter("list", "param", "", null);
        String[] titles = {"all", "your", "base"};
        String[] values = {"a", "y", "b"};
        param.setList(titles, values);
        param.setAttribute("maximum", "2");
        Assert.assertTrue(param.validate("y"));
        Assert.assertTrue(param.validate("b,y"));
        Assert.assertFalse(param.validate("b,x"));
        Assert.assertFalse(param.validate("b,y,a"));
    }
    
    public static void main(String[] args) throws Exception {
        ACDDescriptionTest test = new ACDDescriptionTest();
        test.setUp();
        test.testACDDescription();
        test.testACDParameter();
    }
}