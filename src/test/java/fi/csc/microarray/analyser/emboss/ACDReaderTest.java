package fi.csc.microarray.analyser.emboss;

import org.testng.Assert;

import fi.csc.microarray.analyser.emboss.ACDReader;
import fi.csc.microarray.description.ParsedVVSADL;

public class ACDReaderTest {
    public void run() {
        // FIXME
        ACDReader reader = new ACDReader("water.acd");
        ParsedVVSADL internal = reader.analyseAcd();
        
        Assert.assertEquals(internal.parameters().size(), 2);
        Assert.assertEquals(internal.inputs().size(), 3);
        Assert.assertEquals(internal.outputs().size(), 1);
    }
    
    public static void main(String[] args) throws Exception {
        new ACDReaderTest().run();
    }
}

