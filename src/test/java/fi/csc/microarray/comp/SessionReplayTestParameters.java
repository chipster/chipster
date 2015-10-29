package fi.csc.microarray.comp;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

public class SessionReplayTestParameters {


	@Parameter(names = {"-help", "--help", "-h"}, help = true)
	private boolean help;



    @Parameters(commandDescription = "Run test sessions")
    public static class CommandRun {
        @Parameter(names = "-username", description = "Username", required = true)
        public String username;

        @Parameter(names = "-password", description = "Password", required = true)
        public String password;

        @Parameter(names = "-config", description = "Configuration url", required = true)
        public String config;

        @Parameter(names = "-sessions", description = "Sessions dir", required = true)
        public String sessions;

        @Parameter(names = "-output", description = "Output dir")
        public String output;
    	
    }

    @Parameters(commandDescription = "Generate from previous results")
    public static class CommandGenerate {
        @Parameter(names = "-generate-from", description = "Generate summary files from previous results dir", required = true)
        public String generateFrom;
    	
        @Parameter(names = "-output", description = "Output dir")
        public String output;

    }

    
    
}

