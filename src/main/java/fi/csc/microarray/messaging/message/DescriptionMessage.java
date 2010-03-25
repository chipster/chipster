package fi.csc.microarray.messaging.message;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import fi.csc.microarray.util.XmlUtil;

/**
 * Message for sending module information and tool
 * descriptions.
 * 
 * @author naktinis
 *
 */
public class DescriptionMessage extends ChipsterMessage {
    
    private final static String KEY_MODULE = "module";
    
    private Document moduleXml;
    private String moduleName;
    private Element moduleElement;
    private List<Category> categories = new LinkedList<Category>();
    
    public DescriptionMessage(String moduleName) throws ParserConfigurationException {
        setModuleName(moduleName);
        
        // Start constructing the XML
        moduleXml = XmlUtil.newDocument();
        moduleElement = moduleXml.createElement("module");
        moduleElement.setAttribute("name", moduleName);
        moduleXml.appendChild(moduleElement);
    }
    
    private void setModuleName(String moduleName) {
        this.moduleName = moduleName;
    }
    
    public String getModuleName() {
        return moduleName;
    }
    
    public void addConfString(String configuration) {
        // TODO ADD to XML
    }
    
    public void addCategory(Category category) {
        Element categoryElement = moduleXml.createElement("category");
        categoryElement.setAttribute("name", category.getName());
        String colorString = Integer.toHexString(category.getColor().getRGB());
        colorString = "#" + colorString.substring(2, colorString.length());
        categoryElement.setAttribute("color", colorString);
        
        for (Tool tool : category.getTools()) {
            Element toolElement = moduleXml.createElement("tool");
            toolElement.setAttribute("name", tool.getName());
            toolElement.setTextContent(tool.getDescription());
            categoryElement.appendChild(toolElement);
        }
        
        moduleElement.appendChild(categoryElement);
    }
    
    public List<Category> getCategories() {
        return categories;
    }
    
    @Override
    public void unmarshal(MapMessage from) throws JMSException {
        super.unmarshal(from);
        this.moduleXml = XmlUtil.stringToXML(from.getString(KEY_MODULE));
    }
    
    @Override
    public void marshal(MapMessage to) throws JMSException {
        super.marshal(to);
        to.setString(KEY_MODULE, XmlUtil.xmlToString(moduleXml));
    }
    
    /**
     * Tool category. Contains several related tools.
     */
    public static class Category {
        private String name;
        private List<Tool> tools = new LinkedList<Tool>();
        private Color color;
        
        public Category(String name, String color) {
            this.name = name;
            this.color = Color.decode(color);
        }
        
        public String getName() {
            return name;
        }
        
        public Color getColor() {
            return color;
        }
        
        public void addTool(String name, String description) {
            tools.add(new Tool(name, description));
        }
        
        public List<Tool> getTools() {
            return tools;
        }
    }
    
    /**
     * Single analysis tool.
     */
    public static class Tool {
        private String name;
        private String description;
        
        public Tool(String name, String description) {
            this.name = name;
            this.description = description;
        }
        
        public String getName() {
            return name;
        }
        
        public String getDescription() {
            return description;
        }
    }
}
