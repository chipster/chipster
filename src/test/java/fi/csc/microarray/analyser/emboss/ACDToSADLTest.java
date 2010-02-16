package fi.csc.microarray.analyser.emboss;

import java.io.File;

import org.testng.Assert;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.analyser.emboss.ACDToSADL;

public class ACDToSADLTest {
    String path = "src/test/resources/";
    
    public void run() {

        ACDToSADL reader = new ACDToSADL(new File(path + "water.acd"));
        SADLDescription internal = reader.analyseAcd();
        
        Assert.assertEquals(internal.parameters().size(), 2);
        Assert.assertEquals(internal.inputs().size(), 3);
        Assert.assertEquals(internal.outputs().size(), 1);
    }
    
    public static void main(String[] args) throws Exception {
        new ACDToSADLTest().run();
    }
}