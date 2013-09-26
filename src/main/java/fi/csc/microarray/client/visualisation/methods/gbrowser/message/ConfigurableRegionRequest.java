package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Map;

public class ConfigurableRegionRequest extends DataRequest {

	private Map<String, Object> configuration;

	public ConfigurableRegionRequest(Region region, Map<String, Object> configuration) {
	
		super(region, null, new DataStatus());
		
		this.configuration = configuration;
	}
	
	public Object getConfiguration(String key) {
		return configuration.get(key);
	}
}
