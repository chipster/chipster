package fi.csc.microarray.client.cli;

import java.io.File;
import java.io.IOException;

import org.joda.time.DateTime;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.auth.ClientLoginListener;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;

public class CliClientApplication extends ClientApplication {

	private boolean verbose;
	private boolean quiet;
	
	private DateTime lastLogin;

	public CliClientApplication(SimpleAuthenticationRequestListener auth,
			boolean verbose, boolean quiet) {
		super(auth);
		this.verbose = verbose;
		this.quiet = quiet;		
	}

	@Override
	public void reportExceptionThreadSafely(Exception e) {
		e.printStackTrace();
		System.exit(1);
	}

	@Override
	public void reportException(Exception e) {
		e.printStackTrace();
		System.exit(1);
	}

	@Override
	public void reportTaskError(Task job) throws MicroarrayException {
		System.err.println("\nTask error: " + job.getErrorMessage());
		// no need to exit, CliClient will continue when the job is removed from TaskExecutor
	}
	
	long lastInitMsgTime = System.currentTimeMillis();

	@Override
	public void reportInitialisationThreadSafely(String report, boolean newline) {
		if (verbose) {
			if (newline) {
				System.out.println(report);
			} else {
				System.out.print(report);				
			}
			
//			System.out.println("\t" + (System.currentTimeMillis() - lastInitMsgTime) + "ms\t");
//			lastInitMsgTime = System.currentTimeMillis();
		}
	}

	@Override
	public void showDialog(String title, String message, String details,
			Severity severity, boolean modal) {
		System.out.println(title);
		System.out.println(message);
		System.out.println(details);
		System.exit(1);
	}

	@Override
	public void showDialog(String title, String message, String details,
			Severity severity, boolean modal,
			DetailsVisibility detailsVisibility, PluginButton button) {
		showDialog(title, message, details, severity, modal);
	}

	@Override
	public void showDialog(String title, String message, String details,
			Severity severity, boolean modal,
			DetailsVisibility detailsVisibility, PluginButton button,
			boolean feedBackEnabled) {
		showDialog(title, message, details, severity, modal);
	}

	@Override
	public void runBlockingTask(String taskName, Runnable runnable) {
		if (!quiet) {
			System.out.println("Please wait while " + taskName);
		}
		new Thread(runnable).start(); 
	}

	@Override
	public void initialiseGUIThreadSafely(File mostRecentDeadTempDirectory)
			throws MicroarrayException, IOException {
	}
	
	/* Very complicated way of just closing the program when authentication fails. 
	 * Otherwise server keeps asking for authentication and SimpleAuth keeps sending
	 * wrong credentials.
	 * 
	 * (non-Javadoc)
	 * @see fi.csc.microarray.client.ClientApplication#getAuthenticationRequestListener()
	 */
	@Override
	protected AuthenticationRequestListener getAuthenticationRequestListener() {

		final AuthenticationRequestListener authenticator = super.getAuthenticationRequestListener();

		authenticator.setLoginListener(new ClientLoginListener() {

			public void firstLogin() {
				// assume that authentication failed if it is less than 10 seconds since last authentication request 
				if (lastLogin != null && lastLogin.isAfter(DateTime.now().minus(10000))) {
					System.err.println("login failed");
					System.exit(1);
				}
				
				lastLogin = DateTime.now();				
			}

			public void loginCancelled() {
				System.exit(1);
			}
		});

		return authenticator;
	}
}
