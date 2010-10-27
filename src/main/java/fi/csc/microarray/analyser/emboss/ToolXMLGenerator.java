package fi.csc.microarray.analyser.emboss;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

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
    private List<String> bottomGroups = new LinkedList<String>();
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
        ignoredGroups.add("Ontology:EDAM");

        // less important groups which should be last when sorting groups
        bottomGroups.add("Enzyme Kinetics");
        bottomGroups.add("Feature tables");
        
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
        colors.put("alignment:multiple", "#d59f45");
        colors.put("display", "#e7881c");
        colors.put("edit", "#d53833");
        colors.put("hmm", "#80a3b7");
        colors.put("information", "#0177b7");
        colors.put("nucleic", "#629a9b");
        colors.put("phylogeny", "#a49900");
        colors.put("protein", "#83010b");
        colors.put("enzyme kinetics", "#c0d2de");
        colors.put("feature tables", "#c0d2de");
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
        
        // sort groups
        LinkedList<String> sortedGroups = new LinkedList<String>(groupsMap.keySet());
        Collections.sort(sortedGroups);
        
        // drop less important categories to the bottom
        for (String group : bottomGroups) {
        	if (sortedGroups.contains(group)) {
        		sortedGroups.remove(group);
        		sortedGroups.add(group);
        	}
        }
        
        // generate xml
        for (String group : sortedGroups) {
        	// Create category element
            Element category = doc.createElement("category");
            category.setAttribute("name", group.substring(0,1).toUpperCase() + group.substring(1));
            String colorKey = group.trim().toLowerCase();
            if (colors.containsKey(colorKey)) {
            	category.setAttribute("color", colors.get(colorKey));
            } else {
            	colorKey =  group.split(":")[0].trim().toLowerCase();
                if (colors.containsKey(colorKey)) {
                    category.setAttribute("color", colors.get(colorKey));
                }
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
    	if (args.length != 2) {
    		System.out.println("Please provide <directory for acd files> and <target xml file> as arguments.");
    		return;
    	}
    	
    	new ToolXMLGenerator(args[0], args[1]).generate();
    }
}
