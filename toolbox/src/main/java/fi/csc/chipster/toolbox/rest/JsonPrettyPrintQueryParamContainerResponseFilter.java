package fi.csc.chipster.toolbox.rest;

import java.io.IOException;

import javax.ws.rs.container.ContainerRequestContext;
import javax.ws.rs.container.ContainerResponseContext;
import javax.ws.rs.core.MultivaluedMap;
import javax.ws.rs.ext.Provider;
import com.fasterxml.jackson.jaxrs.cfg.ObjectWriterInjector;


@Provider
public class JsonPrettyPrintQueryParamContainerResponseFilter implements javax.ws.rs.container.ContainerResponseFilter {

    private static final String QUERY_PARAM_PRETTY = "pretty";

    @Override
    public void filter(
            ContainerRequestContext requestContext,
            ContainerResponseContext responseContext) throws IOException {

        MultivaluedMap<String, String> queryParams = requestContext.getUriInfo().getQueryParameters();

        if (queryParams.containsKey(QUERY_PARAM_PRETTY)) {
        	ObjectWriterInjector.set(new IndentingModifier());
        }
    }
}
