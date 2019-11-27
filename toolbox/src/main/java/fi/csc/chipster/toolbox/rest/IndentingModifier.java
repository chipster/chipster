package fi.csc.chipster.toolbox.rest;

import java.io.IOException;

import javax.ws.rs.core.MultivaluedMap;


import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.jaxrs.cfg.EndpointConfigBase;
import com.fasterxml.jackson.jaxrs.cfg.ObjectWriterModifier;

public class IndentingModifier extends ObjectWriterModifier {

    @Override
    public ObjectWriter modify(
            EndpointConfigBase<?> endpoint,
            MultivaluedMap<String, Object> responseHeaders,
            Object valueToWrite,
            ObjectWriter w,
            JsonGenerator g) throws IOException {

    	g.useDefaultPrettyPrinter();

        return w;
    }
}