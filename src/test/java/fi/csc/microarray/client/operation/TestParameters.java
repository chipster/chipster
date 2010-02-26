package fi.csc.microarray.client.operation;

import java.io.IOException;

import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.client.operation.parameter.EnumParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.operation.parameter.ParameterPanel;
import fi.csc.microarray.client.operation.parameter.ToolParameterPanel;
import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class TestParameters {
    
    @BeforeSuite
    protected static void setUp() throws Exception {
        DirectoryLayout.initialiseClientLayout();
    }
    
    @Test
    public static void runTest() {
        OperationCategory category = new OperationCategory("Testational");
        OperationDefinition definition = new OperationDefinition("Testation",
            category, "Testationing", false);
        Operation operation;
        try {           
            // Prepare the parameters
            SelectionOption[] options = new SelectionOption[2];
            options[0] = new SelectionOption("I'm worth clicking", "i");
            options[1] = new SelectionOption("Click me", "m");
            Parameter param = new EnumParameter("enumpar", "Testy parameter", options, 0, 1, 2);
            definition.addParameter(param);
            
            // Initialize some mock context
            DataBean[] dataBeans = new DataBean[0];
            operation = new Operation(definition, dataBeans);
            
            // Create the panel which maps Parameters to InputComponents
            ToolParameterPanel panel = new ToolParameterPanel(operation, null);
        } catch (MicroarrayException e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        try {
            setUp();
            runTest();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
