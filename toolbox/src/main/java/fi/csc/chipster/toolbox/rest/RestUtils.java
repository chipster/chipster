package fi.csc.chipster.toolbox.rest;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.glassfish.grizzly.GrizzlyFuture;
import org.glassfish.grizzly.http.server.HttpServer;
import org.glassfish.jersey.CommonProperties;
import org.glassfish.jersey.server.ResourceConfig;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.datatype.jsr310.JavaTimeModule;
import com.fasterxml.jackson.jaxrs.json.JacksonJaxbJsonProvider;

public class RestUtils {
	private static Logger logger = LogManager.getLogger();

	private static ObjectMapper objectMapperDefault;

	public static ObjectMapper getObjectMapper() {
		if (objectMapperDefault == null) {

			// separate instance, because configuration may not be thread safe
			objectMapperDefault = getNewObjectMapper();
		}
		return objectMapperDefault;
	}

	public static ObjectMapper getNewObjectMapper() {
		return new ObjectMapper().registerModule(new JavaTimeModule())
				.configure(SerializationFeature.WRITE_DATES_AS_TIMESTAMPS, false);
	}

	public static ResourceConfig getResourceConfig() {

		ResourceConfig rc = new ResourceConfig()
				/*
				 * Disable auto discovery so that we can decide what we want to register and
				 * what not. Don't register JacksonFeature, because it will register
				 * JacksonMappingExceptionMapper, which annoyingly swallows response's
				 * JsonMappingExceptions. Register directly the JacksonJaxbJsonProvider which is
				 * enough for the actual JSON conversion (see the code of JacksonFeature).
				 */
				.property(CommonProperties.FEATURE_AUTO_DISCOVERY_DISABLE, true).register(JacksonJaxbJsonProvider.class)
//				.register(JavaTimeObjectMapperProvider.class)
//				// register all exception mappers
//				.packages(NotFoundExceptionMapper.class.getPackage().getName())
//				// enable the RolesAllowed annotation
//				.register(RolesAllowedDynamicFeature.class)
				.register(JsonPrettyPrintQueryParamContainerResponseFilter.class);

		return rc;
	}

	public static void shutdown(String name, HttpServer httpServer) {

		if (httpServer == null) {
			logger.warn("can't shutdown " + name + ", the server is null");
			return;
		}
		GrizzlyFuture<HttpServer> future = httpServer.shutdown();
		try {
			// wait for server to shutdown, otherwise the next test set will print ugly log
			// messages
			try {
				future.get(3, TimeUnit.SECONDS);
			} catch (TimeoutException e) {
				logger.warn(name + " server didn't stop gracefully");
				httpServer.shutdownNow();
			}
		} catch (InterruptedException | ExecutionException e) {
			logger.warn("failed to shutdown the server " + name, e);
		}
	}

}
