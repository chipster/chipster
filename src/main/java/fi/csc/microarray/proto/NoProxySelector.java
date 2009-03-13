package fi.csc.microarray.proto;


import java.io.IOException;
import java.net.Proxy;
import java.net.ProxySelector;
import java.net.SocketAddress;
import java.net.URI;
import java.util.LinkedList;
import java.util.List;

public class NoProxySelector extends ProxySelector {
	
	@Override
	public void connectFailed(URI uri, SocketAddress sa, IOException ioe) {
		// we are not interested in this
	}

	@Override
	public List<Proxy> select(URI uri) {
		LinkedList<Proxy> proxies = new LinkedList<Proxy>();
		proxies.add(Proxy.NO_PROXY);
		return proxies;
	}

}
