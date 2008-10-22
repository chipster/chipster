package fi.csc.microarray.proto;

import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;

import fi.csc.microarray.MicroarrayMain;

/**
 * Usage on corona:
 * ./jvsc/jdk1.3.1_solaris8_sparc -pidfile nami.pid -debug -cp microarray-current-development.jar fi.csc.microarray.proto.NamiDaemon analyser
 * 
 * @author akallio
 */
public class NamiDaemon implements Daemon {

	String[] args = null;
	
	public void init(DaemonContext context) throws Exception {	
		this.args = context.getArguments();
	}

	public void start() throws Exception {
		MicroarrayMain.main(args);		
	}

	public void stop() throws Exception {		
	}

	public void destroy() {		
	}
}
