package fi.csc.chipster.web.adminweb;

import java.io.File;

import org.eclipse.jetty.security.ConstraintMapping;
import org.eclipse.jetty.security.ConstraintSecurityHandler;
import org.eclipse.jetty.security.HashLoginService;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.server.handler.DefaultHandler;
import org.eclipse.jetty.server.handler.HandlerCollection;
import org.eclipse.jetty.util.security.Constraint;
import org.eclipse.jetty.util.security.Password;
import org.eclipse.jetty.webapp.WebAppContext;

public class StandaloneAdminWeb {
	public static void main(String args[]) throws Exception {
		org.eclipse.jetty.server.Server adminServer = new org.eclipse.jetty.server.Server();
		ServerConnector connector = new ServerConnector(adminServer);
		connector.setPort(8083);
		adminServer.setConnectors(new Connector[]{ connector });
		
		Constraint constraint = new Constraint();
		constraint.setName(Constraint.__BASIC_AUTH);
		constraint.setRoles(new String[] {"admin_role"});
		constraint.setAuthenticate(true);
		
		ConstraintMapping cm = new ConstraintMapping();
		cm.setConstraint(constraint);
		cm.setPathSpec("/*");
		
		HashLoginService loginService = new HashLoginService("Please enter Chipster Admin username and password");
		loginService.update("chipster", 
				new Password("chipster"), 
				new String[] {"admin_role"});
		
		ConstraintSecurityHandler sh = new ConstraintSecurityHandler();
		sh.setLoginService(loginService);
		sh.addConstraintMapping(cm);
		
		WebAppContext context = new WebAppContext();
		File war = new File("../chipster/dist/admin-web.war");
		//File war = new File("webapps/admin-web.war");
		context.setWar(war.getAbsolutePath());
		System.out.println(war.getAbsolutePath());
        context.setContextPath("/");
				
        context.setHandler(sh);
		HandlerCollection handlers = new HandlerCollection();
		handlers.setHandlers(new Handler[] {context, new DefaultHandler()});
				
		adminServer.setHandler(handlers);
        adminServer.start();
	}
}
