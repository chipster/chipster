package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;

import fi.csc.microarray.analyser.emboss.ACDDescription;

public class ACDDescriptionTest {
    String path = "src/test/resources/";
    
    public void run() {
        ACDDescription acd = new ACDDescription(new File(path + "emma.acd"));
        
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
        Assert.assertTrue(!acd.getParameter("maxdiv").isRequired());
        Assert.assertTrue(acd.getParameter("maxdiv").isAdditional());
        Assert.assertTrue(acd.getParameter("sequence").isRequired());
        Assert.assertTrue(!acd.getParameter("sequence").isAdditional());
        Assert.assertTrue(!acd.getParameter("sequence").isAdvanced());
        
        // Test parameter validation
        Assert.assertTrue(acd.getParameter("pwgapopen").validate("0.1"));
        Assert.assertTrue(!acd.getParameter("pwgapopen").validate("-0.1"));
    }
    
    public static void main(String[] args) throws Exception {
        new ACDDescriptionTest().run();
    }
}