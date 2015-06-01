package fi.csc.microarray.util;

import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedHashMap;

import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSocketFactory;
import javax.net.ssl.TrustManager;

/**
 * Create a new SSLSocketFactory for each thread keep a few latest factories in
 * cache.
 * 
 * Sharing a socket factory between threads mostly seems to work fine, but it
 * doesn't promise to be thread safe. When making hundreds of requests to read
 * files partially (in type tagging), the Jetty started to complain weird
 * "Idle timeout expired" messages. Creating a new socket factory every time
 * fixes these symptoms, but affects performance. Probably it breaks connection's
 * keep-alive and forces new SSL handshakes.
 * 
 * Creating a new instance for each thread solves the concurrency problems and
 * the cache fixes the performance implications.
 * 
 * @author klemela
 *
 */
public class SocketFactoryCache {

	/*
	 * By default there are only 5 concurrent HTTP connections in Java
	 * HttpURLConnection, but the genome browser uses many more threads.
	 */
	private static final int CACHE_SIZE = 100;
	private TrustManager[] trustManagers;
	private LinkedHashMap<Thread, SSLSocketFactory> cache = new LinkedHashMap<>();
	private String protocol;

	public SocketFactoryCache(TrustManager[] trustManagers, String protocol) {
		this.trustManagers = trustManagers;
		this.protocol = protocol;
	}

	public SSLSocketFactory getSocketFactoryForThisThread()
			throws NoSuchAlgorithmException, KeyManagementException {
		synchronized (this) {
			if (!cache.containsKey(Thread.currentThread())) {
				// make space
				while (cache.size() > CACHE_SIZE - 1) {
					removeLast(cache);
				}

				SSLContext ctx = SSLContext.getInstance(protocol);
				ctx.init(null,trustManagers, null);
				cache.put(Thread.currentThread(), ctx.getSocketFactory());
			}
			return cache.get(Thread.currentThread());
		}
	}

	private void removeLast(LinkedHashMap<?, ?> map) {
		map.remove(map.keySet().toArray()[map.size() - 1]);
	}
}
