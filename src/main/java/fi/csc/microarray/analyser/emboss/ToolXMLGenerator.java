package fi.csc.microarray.analyser.emboss;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
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
 * <p>
 * TODO: check if a binary file exists for a given
 * acd file.
 * 
 * @author naktinis
 *
 */
public class ToolXMLGenerator {

    // There are problems with 3 acd files: intconv.acd, complex.acd
    // and ensembltest.acd.
    // The first two have 'relations: ' instead of the correct
    // 'relations: ""'. The third one (ensembltest.acd)
    // had 'default: 25' instead of 'default: "25"'
    // NOTE: they have been reported and probably fixed in Emboss' cvs

    private File acdDir;
    private File outFile;
    private Document doc = null;
    private static final HashSet<String> ignoredGroups = new HashSet<String>();
    private static final HashSet<String> ignoredPrograms = new HashSet<String>();
    private LinkedHashMap<String, LinkedList<String>> groupsMap = new LinkedHashMap<String, LinkedList<String>>();
    private HashMap<String, String> colors = new HashMap<String, String>();

    public ToolXMLGenerator(String acdDirPath, String outputFilePath) {
        acdDir = new File(acdDirPath);
        outFile = new File(outputFilePath);

        // Top-level groups we are not interested in
        // See:
        // http://emboss.sourceforge.net/developers/acd/syntax.html#sect2214
        ignoredGroups.add("");
        ignoredGroups.add("Acd");
        ignoredGroups.add("Menus");
        ignoredGroups.add("Test");
        ignoredGroups.add("Utils");

        // Separate programs that are broken or we are not
        // interested in (they might be in several groups)
        
        // Programs that are not present
        ignoredPrograms.add("finddb");
        ignoredPrograms.add("showdball");
        ignoredPrograms.add("seqxrefall");
        ignoredPrograms.add("showid");
        ignoredPrograms.add("idtell");
        ignoredPrograms.add("idxref");
        ignoredPrograms.add("dbquery");
        ignoredPrograms.add("prima");
        ignoredPrograms.add("primers");
        ignoredPrograms.add("newcoils");
        ignoredPrograms.add("ememe");
        
        // Databases for these programs not supported
        ignoredPrograms.add("tfscan");
        
        // FIXME Wait until these descriptions fixed in the new EMBOSS version
        ignoredPrograms.add("dbxstat");
        ignoredPrograms.add("emiraest");
        ignoredPrograms.add("ememetext");
        ignoredPrograms.add("emira");
        ignoredPrograms.add("vrnafoldpf");        
        
        colors.put("alignment", "#e7df70");
        colors.put("display", "#d59f45");
        colors.put("edit", "#e7881c");
        colors.put("enzyme kinetics", "#83010b");
        colors.put("feature tables", "#83010b");
        colors.put("information", "#80a3b7");
        colors.put("nucleic", "#0177b7");
        colors.put("phylogeny", "#629a9b");
        colors.put("protein", "#a49900");
        colors.put("hmm", "#d53833");
    }

    /**
     * Generate an XML tree storing information about the tools available in
     * EMBOSS package.
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
        Element module = doc.createElement("module");
        module.setAttribute("name", "sequence");
        doc.appendChild(module);
        
        for (File acdFile : acdDir.listFiles()) {
            System.out.println(acdFile.getName());
            try {
                if (acdFile.getName().endsWith(".acd")) {
                    // Feed file content to the parser
                    BufferedInputStream inputStream =
                        new BufferedInputStream(new FileInputStream(acdFile));
                    final byte [] bytes = new byte[(int) acdFile.length()];
                    inputStream.read(bytes);
                    
                    // Check if this program is not ignored
                    if (ignoredPrograms.contains(
                            acdFile.getName().substring(0, acdFile.getName().indexOf(".")))) {
                        continue;
                    }
                    
                    ACDDescription acd = new ACDDescription(acdFile);
                    
                    // Check if there are any interesting groups for this application
                    HashSet<String> acdGroups =
                        new HashSet<String>(getTopLevelGroups(acd.getGroups()));
                    acdGroups.removeAll(ignoredGroups);
                    LinkedList<String> acdGroupList = new LinkedList<String>();
                    acdGroupList.addAll(acdGroups);
                    
                    // Check if we are interested in this group
                    if (acdGroups.size() > 0) {
                        // Add this application to group map so we can sort it 
                        // and add to xml later.
                        // TODO: add all groups, not only the first one
                        if (!groupsMap.containsKey(acd.getGroups().get(0))) {
                            groupsMap.put(acd.getGroups().get(0), new LinkedList<String>());
                        }
                        groupsMap.get(acd.getGroups().get(0)).add(acdFile.getName());
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        // Sort groups and generate XML
        LinkedList<String> sortedGroups = new LinkedList<String>(groupsMap.keySet());
        Collections.sort(sortedGroups);
        for (String group : sortedGroups) {
            // Create category element
            Element category = doc.createElement("category");
            category.setAttribute("name", group.substring(0,1).toUpperCase() + group.substring(1));
            String groupNormal = group.split(":")[0].trim().toLowerCase();
            if (colors.containsKey(groupNormal)) {
                category.setAttribute("color", colors.get(groupNormal));
            }
            module.appendChild(category);
            LinkedList<String> sortedApps = groupsMap.get(group);
            Collections.sort(sortedApps);
            for (String appName : sortedApps) {
                
                // Make an entry in the XML tree
                Element tool = doc.createElement("tool");
                Element resource = doc.createElement("resource");
                
                // The filename is more important than application name
                // according to ACD specification
                resource.setTextContent(appName);
                tool.setAttribute("runtime", "EMBOSS");
                tool.appendChild(resource);
                category.appendChild(tool);
            }
        }

        saveToFile(outFile);
    }

    /**
     * Take a list of groups and extract the top-level names.
     * 
     * Before calling this we might have:
     * 
     * <pre>
     * ["Phylogeny:Consensus", "Phylogeny:Misc"]
     * </pre>
     * 
     * After this we have:
     * 
     * <pre>
     * ["Phylogeny", "Phylogeny"]
     * </pre>
     * 
     * @param groupString
     * @return
     */
    private LinkedList<String> getTopLevelGroups(
            LinkedList<String> origGroupList) {
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
            serializer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
            serializer.setOutputProperty(OutputKeys.INDENT, "yes");

            Source source = new DOMSource(doc);
            Result result = new StreamResult(file);

            serializer.transform(source, result);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new ToolXMLGenerator("/home/naktinis/acd",
                "debug-base-dir/conf/sequence-module.xml").generate();
    }
}
