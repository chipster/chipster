package fi.csc.chipster.toolbox;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Collection;
import java.util.List;

import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;

public class ToolboxRestUtils {

	public static String asJson(Object obj) throws JsonGenerationException, JsonMappingException, IOException {
		// using Jackson library
		StringWriter writer = new StringWriter();
		// support for LocalDateTime
		ObjectMapper mapper = new ObjectMapper();
		mapper.writeValue(writer, obj);
		return writer.toString();
	}

	public static <T> T parseJson(Class<T> obj, String json)
			throws JsonParseException, JsonMappingException, IOException {
		return parseJson(obj, json, true);
	}

	public static <T> T parseJson(Class<T> obj, String json, boolean failOnUnknownProperties)
			throws JsonParseException, JsonMappingException, IOException {
		// using Jackson library
		StringReader reader = new StringReader(json);
		ObjectMapper mapper = new ObjectMapper();
		mapper.configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, failOnUnknownProperties);
		return mapper.readValue(reader, obj);
	}

	@SuppressWarnings("rawtypes")
	public static List parseJson(Class<? extends Collection> collectionType, Class<?> itemType, String json)
			throws JsonParseException, JsonMappingException, IOException {
		return parseJson(collectionType, itemType, json, true);
	}

	@SuppressWarnings("rawtypes")
	public static List parseJson(Class<? extends Collection> collectionType, Class<?> itemType, String json,
			boolean failOnUnknownProperties) throws JsonParseException, JsonMappingException, IOException {
		// using Jackson library
		StringReader reader = new StringReader(json);
		// support for LocalDateTime
		ObjectMapper mapper = new ObjectMapper();
		mapper.configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, failOnUnknownProperties);
		return mapper.readValue(reader, mapper.getTypeFactory().constructCollectionType(collectionType, itemType));
	}

}
