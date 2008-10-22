package fi.csc.microarray.util.rest;

import java.io.IOException;

import javax.servlet.Filter;
import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletRequestWrapper;

import org.mortbay.log.Log;

public class FilenameGuardFilter implements Filter {

	public void destroy() {
		// nothing to destroy		
	}

	public void init(FilterConfig config) throws ServletException {
		// nothing to init		
	}

	public void doFilter(ServletRequest request, ServletResponse response, FilterChain chain) throws IOException, ServletException {

		if (request instanceof HttpServletRequest) {
			HttpServletRequest httpRequest = (HttpServletRequest)request;
			GuardedHttpServletRequest guardedRequest = new GuardedHttpServletRequest(httpRequest);
			chain.doFilter(guardedRequest, response);
		} else {
			chain.doFilter(request, response);
		}		
	}
	
	private static class GuardedHttpServletRequest extends HttpServletRequestWrapper {
		
		public GuardedHttpServletRequest(HttpServletRequest httpRequest) {
			super(httpRequest);
		}

		private String guard(String filename) {
			if (Log.isDebugEnabled()) {
				Log.debug("GUARDed " + filename + " to "+ filename.replace(":", "_"));
			}			
			return filename.replace(":", "_");
		}
		
		@Override
		public String getParameter(String name) {
			if (name.equals("Destination")) {
				return guard(super.getParameter(name));
			} else {
				return super.getParameter(name);
			}
		}
		
		@Override
		public String getPathInfo() {
			return guard(super.getPathInfo());
		}
		
		@Override
		public String getPathTranslated() {
			return guard(super.getPathTranslated());
		}
		
		@Override
		public String getRequestURI() {
			return guard(super.getRequestURI());
		}
	}
}
