package fi.csc.microarray.analyser.emboss;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.HashSet;
import java.util.LinkedList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Generate an XML file with information about
 * tools available in EMBOSS package.
 * 
 * @author naktinis
 *
 */
public class ToolXMLGenerator {
    
    // TODO
    // There are problems with 3 acd files: intconv.acd, complex.acd
    // and ensembltest.acd.
    // The first two have 'relations: ' instead of the correct
    // 'relations: ""'. The third one (ensembltest.acd)
    // had 'default: 25' instead of 'default: "25"' 
    
    private File acdDir;
    private File outFile;
    private Document doc = null;
    private static final HashSet<String> ignoredGroups = new HashSet<String>();
    
    public ToolXMLGenerator(String acdDirPath, String outputFilePath) {
        acdDir = new File(acdDirPath);
        outFile = new File(outputFilePath);
        
        // Top-level groups we are not interested in
        // See: http://emboss.sourceforge.net/developers/acd/syntax.html#sect2214
        ignoredGroups.add("");
        ignoredGroups.add("Acd");
        ignoredGroups.add("Menus");
        ignoredGroups.add("Test");
        ignoredGroups.add("Utils");
    }
    
    /**
     * Generate an XML tree storing information about
     * the tools available in EMBOSS package.
     */
    public void generate() {
        try {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            doc = builder.newDocument();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }
        
        // Create a node for all tools
        Element tools = doc.createElement("tools");
        doc.appendChild(tools);
        
        // Read all acd files in given directory
        for (File acdFile : acdDir.listFiles()) {
            try {
                if (acdFile.getName().endsWith(".acd")) {
                    // Feed file content to the parser
                    BufferedInputStream inputStream =
                        new BufferedInputStream(new FileInputStream(acdFile));
                    final byte [] bytes = new byte[(int) acdFile.length()];
                    inputStream.read(bytes);
                    ACDDescription acd = new ACDDescription(acdFile);
                    
                    HashSet<String> acdGroups =
                        new HashSet<String>(getTopLevelGroups(acd.getGroups()));
                    acdGroups.removeAll(ignoredGroups);
                    
                    // Check if we are interested in this group
                    if (acdGroups.size() > 0) {
                        // Make an entry in the XML tree
                        Element tool = doc.createElement("tool");
                        
                        // The filename is more important than application name
                        // according to ACD specification
                        tool.setTextContent(acdFile.getName());
                        tool.setAttribute("runtime", "EMBOSS");
                        tools.appendChild(tool);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        saveToFile(outFile);
    }
    
    /**
     * Take a list of groups and extract the top-level names.
     * 
     * Before calling this we might have:<pre>
     * ["Phylogeny:Consensus", "Phylogeny:Misc"]</pre>
     * 
     * After this we have:<pre>
     * ["Phylogeny", "Phylogeny"]</pre>
     * 
     * @param groupString
     * @return
     */
    private LinkedList<String> getTopLevelGroups(LinkedList<String> origGroupList) {
        LinkedList<String> groupList = new LinkedList<String>();
        for (String group : origGroupList) {
            groupList.add(group.split(":")[0].trim());
        }
        return groupList;
    }
    
    private void saveToFile(File file) {
        // Make XML writable to a file
        TransformerFactory tf = TransformerFactory.newInstance();
        Transformer serializer;
        try {
            serializer = tf.newTransformer();
            serializer.setOutputProperty(OutputKeys.ENCODING,"UTF-8");
            serializer.setOutputProperty(OutputKeys.INDENT,"yes");
            
            Source source = new DOMSource(doc);
            Result result = new StreamResult(file);
            
            serializer.transform(source, result);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        new ToolXMLGenerator("/opt/EMBOSS-6.2.0/emboss/acd",
                             "debug-base-dir/conf/emboss-tools.xml").generate();
    }
}
