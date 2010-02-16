package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;

import fi.csc.microarray.analyser.emboss.ACD;

public class ACDTest {
    String path = "src/test/resources/";
    
    public void run() {
        ACD acd = new ACD();
        acd.fromFile(new File(path + "emma.acd"));
              
        Assert.assertEquals(acd.getParameters("input", null, false).size(), 5);
        Assert.assertEquals(acd.getParameters("input", null, true).size(), 11);
        Assert.assertEquals(acd.getParameters("input", "matrixsection", false).size(), 3);
        Assert.assertEquals(acd.getParameters("additional", null, true).size(), 15);
        
        Assert.assertEquals(acd.getParameter("maxdiv").getAttribute("default"), "30");
        Assert.assertEquals(acd.getParameter("dendoutfile").getAttribute("extension"), "dnd");
    }
    
    public static void main(String[] args) throws Exception {
        new ACDTest().run();
    }
}