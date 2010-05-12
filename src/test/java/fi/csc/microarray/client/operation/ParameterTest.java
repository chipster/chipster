package fi.csc.microarray.client.operation;

import java.awt.BorderLayout;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.client.operation.parameter.EnumParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.operation.parameter.ToolParameterPanel;
import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Test client-side parameter creation.
 * 
 * @author naktinis
 */
public class ParameterTest {
    
    private ToolParameterPanel panel;
    private Parameter paramMulti;
    private Parameter paramSingle;
    
    @BeforeSuite
    protected void setUp() throws Exception {
        DirectoryLayout.initialiseClientLayout();
    }
    
    @BeforeClass
    public void prepareComponents() {
        OperationCategory category = new OperationCategory("Testational");
        OperationDefinition definition = new OperationDefinition("Testation id", null,
            category, "Testationing", false);
        Operation operation;
        try {
            Name[] options = {Name.createName("a"), Name.createName("b")};
            String[] defaults = {"a"};
            Parameter p = Parameter.createInstance(Name.createName("list"), ParameterType.ENUM,
                                                   options, "This is list", "1", "1", defaults, false);
            Assert.assertFalse(p.isOptional());
            
            // Prepare the multi-select parameter
            SelectionOption[] optionsMulti = new SelectionOption[3];
            optionsMulti[0] = new SelectionOption("I'm worth clicking", "i");
            optionsMulti[1] = new SelectionOption("Click me", "m");
            optionsMulti[2] = new SelectionOption("Or me", "o");
            List<SelectionOption> defaultOptions = new LinkedList<SelectionOption>();
            defaultOptions.add(optionsMulti[0]);
            defaultOptions.add(optionsMulti[2]);
            paramMulti = new EnumParameter("multi", "Enum parameter", "Enum parameter", optionsMulti, defaultOptions, 1, 2);
            definition.addParameter(paramMulti);
            
            // Prepare a single-select parameter
            SelectionOption[] optionsSingle = new SelectionOption[2];
            optionsSingle[0] = new SelectionOption("I'm worth choosing", "i");
            optionsSingle[1] = new SelectionOption("Choose me", "m");
            defaultOptions = new LinkedList<SelectionOption>();
            defaultOptions.add(optionsSingle[1]);
            paramSingle = new EnumParameter("single", "Enum parameter", "Enum parameter", optionsSingle, defaultOptions, 1, 1);
            definition.addParameter(paramSingle);
            
            // Initialize some mock context
            DataBean[] dataBeans = new DataBean[0];
            operation = new Operation(definition, dataBeans);
            
            // Create the panel which maps Parameters to InputComponents
            panel = new ToolParameterPanel(operation, null);
        } catch (MicroarrayException e) {
            e.printStackTrace();
        }
    }
    
    @Test
    public void runTest() {

        // Check if component created successfully
        Assert.assertEquals(panel.getComponentCount(), 1);         
        
        // Check if default values are correct
        Assert.assertEquals(paramMulti.getValue(), "i,o");
        Assert.assertEquals(paramSingle.getValue(), "m");

    }
    
    /**
     * Display a window containing several UI components
     * that represent parameters. This might need commenting
     * out some lines in actionPerformed methods or having
     * more mock objects.
     */
    private void testVisually() {
        // Create a frame
        JFrame frame = new JFrame("Test frame");

        // Create components and put them in the frame
        frame.getContentPane().add(panel.getComponent(0), BorderLayout.CENTER);

        // Size it and show it
        frame.pack();
        frame.setVisible(true);
    }
    
    public static void main(String[] args) {
        ParameterTest test = new ParameterTest();
        try {
            test.setUp();            
            test.prepareComponents();
            test.runTest();
            test.testVisually();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
