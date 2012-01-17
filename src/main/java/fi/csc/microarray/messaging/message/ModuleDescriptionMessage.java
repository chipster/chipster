package fi.csc.microarray.messaging.message;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.microarray.util.XmlUtil;

/**
 * Message for sending module information and tool
 * descriptions.
 * 
 * @author naktinis
 *
 */
public class ModuleDescriptionMessage extends ChipsterMessage {
    
    private final static String KEY_MODULE = "module";
    private final static String KEY_MODULE_NAME = "module-name";
    
    public Document moduleXml;
    private String moduleName;
    private List<Category> categories = new LinkedList<Category>();
    
    /**
     * Empty constructor (needed for MessageListenerWrap.onMessage)
     */
    public ModuleDescriptionMessage() { }
    
    public ModuleDescriptionMessage(String moduleName) {
    	try {
    		// Start constructing the XML
    		moduleXml = XmlUtil.newDocument();
    		moduleXml.appendChild(moduleXml.createElement("module"));

    		// Set module's name
    		setModuleName(moduleName);
    		getFirstModule().setAttribute("name", moduleName);

    	} catch (ParserConfigurationException e) {
    		throw new RuntimeException(e); // should never happen
    	}
    }
    
    private void setModuleName(String moduleName) {
        this.moduleName = moduleName;
    }
    
    public String getModuleName() {
        return this.moduleName;
    }
    
    public void addConfString(String configuration) {
        // TODO ADD to XML
    }
    
    private Element getFirstModule() {
        return (Element)moduleXml.getElementsByTagName("module").item(0);
    }
    
    /**
     * Store a Category object in the inner XML.
     * 
     * @param category
     */
    public void addCategory(Category category) {
        Element categoryElement = moduleXml.createElement("category");
        categoryElement.setAttribute("name", category.getName());
        String colorString = Integer.toHexString(category.getColor().getRGB());
        colorString = "#" + colorString.substring(2, colorString.length());
        categoryElement.setAttribute("color", colorString);
        categoryElement.setAttribute("hidden", category.isHidden().toString());
        
        for (Tool tool : category.getTools()) {
            Element toolElement = moduleXml.createElement("tool");
            toolElement.setAttribute("helpURL", tool.getHelpURL());
            toolElement.setTextContent(tool.getDescription());
            categoryElement.appendChild(toolElement);
        }
        
        getFirstModule().appendChild(categoryElement);
    }
    
    /**
     * @return all categories contained in the inner XML.
     */
    public List<Category> getCategories() {
        NodeList categoryList = getFirstModule().getElementsByTagName("category");
        for (int i=0; i<categoryList.getLength(); i++) {
            Element categoryElement = (Element) categoryList.item(i);
            Category category = new Category(categoryElement.getAttribute("name"),
                    categoryElement.getAttribute("color"),
                    Boolean.valueOf(categoryElement.getAttribute("hidden")));
            NodeList toolList = categoryElement.getElementsByTagName("tool");
            for (int j=0; j<toolList.getLength(); j++) {
                Element toolElement = (Element) toolList.item(j);
                category.addTool(toolElement.getTextContent(),
                        toolElement.getAttribute("helpURL"));
            }
            categories.add(category);
        }
        return categories;
    }
    
    @Override
    public void unmarshal(MapMessage from) throws JMSException {
        super.unmarshal(from);
        this.setModuleName(from.getStringProperty(KEY_MODULE_NAME));
        this.moduleXml = XmlUtil.stringToXML(from.getString(KEY_MODULE));
    }
    
    @Override
    public void marshal(MapMessage to) throws JMSException {
        super.marshal(to);
        to.setStringProperty(KEY_MODULE_NAME, this.getModuleName());
        to.setString(KEY_MODULE, XmlUtil.xmlToString(moduleXml));
    }
    
    /**
     * Tool category. Contains several related tools.
     */
    public static class Category {
        private String name;
        private List<Tool> tools = new LinkedList<Tool>();
        private Color color;
        private Boolean hidden;
        
        /**
         * @param name - name for this category.
         * @param color - String representing hexidecimal RGB value (e.g. "FF1122")
         */
        public Category(String name, String color, Boolean hidden) {
            this.name = name;
            this.color = Color.decode(color);
            this.hidden = hidden;
        }
        
        public String getName() {
            return name;
        }
        
        public Color getColor() {
            return color;
        }
        
        public Boolean isHidden() {
            return hidden;
        }
        
        public void addTool(String description, String helpURL) {
            tools.add(new Tool(description, helpURL));
        }
        
        public List<Tool> getTools() {
            return tools;
        }
    }
    
    /**
     * Single analysis tool.
     */
    public static class Tool {
        private String description;
        private String helpURL;
        
        public Tool(String description, String helpURL) {
            this.description = description;
            this.helpURL = helpURL;
        }
        
        public String getHelpURL() {
            return helpURL;
        }
        
        public String getDescription() {
            return description;
        }
    }
}
