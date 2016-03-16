package fi.csc.chipster.toolbox;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;

import fi.csc.microarray.description.SADLDescription;

public class ToolboxClientImpl {

    public static void main(String args[]) {
        Client client = ClientBuilder.newClient();
//        System.out.println(config.getString("toolbox"));
        WebTarget webTarget = client.target("http://127.0.0.1:8008/toolbox").path("tools/norm-affy.R/sadl");

//        ToolboxTool tool = webTarget.request(MediaType.APPLICATION_JSON).get(ToolboxTool.class);
        String sadl = webTarget.request().get(String.class);
        System.out.println(sadl);
    }
    

	private String baseUri;
	private Client client;
	
	public ToolboxClientImpl(String toolboxUri) {
		this.baseUri = toolboxUri;
		this.client = ClientBuilder.newClient();;
	}

	public ToolboxTool getTool(String toolId) throws JsonParseException, JsonMappingException, IOException {
		
		WebTarget serviceTarget = client.target(baseUri).path("tools/" + toolId);

		String json = serviceTarget.request(MediaType.APPLICATION_JSON).get(String.class);
		
		ToolboxTool tool = ToolboxRestUtils.parseJson(ToolboxTool.class, json, false);
		
		return tool;
	}

	public HashMap<String, SADLDescription> getTools() throws JsonParseException, JsonMappingException, IOException {
		WebTarget serviceTarget = client.target(baseUri).path("tools");

		String toolsJson = serviceTarget.request(MediaType.APPLICATION_JSON).get(String.class);
		
		@SuppressWarnings("unchecked")
		List<SADLDescription> tools = ToolboxRestUtils.parseJson(List.class, SADLDescription.class, toolsJson, false);
		
		HashMap<String, SADLDescription> map = new HashMap<>();
		
		for (SADLDescription tool : tools) {
			map.put(tool.getName().getID(), tool);
		}

		return map;
	}
	
	public void close() {
		client.close();
	}

}
