package fi.csc.chipster.toolbox;

import java.io.IOException;

import javax.ws.rs.NotFoundException;
import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;

public class ToolboxClientComp {

	private String baseUri;
	private Client client;

	public ToolboxClientComp(String toolboxUri) {
		this.baseUri = toolboxUri;
		this.client = ClientBuilder.newClient();
	}

	public ToolboxTool getTool(String toolId) throws IOException {

		WebTarget serviceTarget = client.target(baseUri).path("tools/" + toolId);

		String json; 
		try {
			json = serviceTarget.request(MediaType.APPLICATION_JSON).get(String.class);
		} catch (NotFoundException nfe) {
			return null;
		}

		ToolboxTool tool = ToolboxRestUtils.parseJson(ToolboxTool.class, json, false);

		return tool;
	}

	public void close() {
		client.close();
	}

	public static void main(String args[]) throws JsonParseException, JsonMappingException, IOException {
		
		ToolboxClientComp toolboxClient = new ToolboxClientComp("http://localhost:8008/toolbox");
		try {
			ToolboxTool tool = toolboxClient.getTool("norm-affy.R");
		
			System.out.println(tool.getSadlString());
		} finally {
			toolboxClient.close();
		}
	}

}
