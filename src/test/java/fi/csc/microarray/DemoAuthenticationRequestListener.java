package fi.csc.microarray;

import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;

public class DemoAuthenticationRequestListener extends SimpleAuthenticationRequestListener {
	public DemoAuthenticationRequestListener() {
		super("testing", "testing");
	}
}
