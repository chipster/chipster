package fi.csc.microarray.comp.java;

import fi.csc.microarray.comp.JobCancelledException;
import fi.csc.microarray.comp.OnDiskCompJobBase;
import fi.csc.microarray.comp.ToolDescription;
import fi.csc.microarray.comp.ToolDescription.ParameterDescription;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;

public abstract class JavaCompJobBase extends OnDiskCompJobBase {

	public static class JavaParameterSecurityPolicy extends ParameterSecurityPolicy {

		private static final int MAX_VALUE_LENGTH = 10000;

		@Override
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {

			// Check parameter size (DOS protection)
			if (value.length() > MAX_VALUE_LENGTH) {
				return false;
			}

			// No need to check content, parameters are passed inside Java Strings
			return true;
		}

		@Override
		public boolean allowUncheckedParameters(ToolDescription toolDescription) {
			return "fi.csc.chipster.tools.common.DownloadFile.java".equals(toolDescription.getID());
		}

	}

	public static JavaParameterSecurityPolicy JAVA_PARAMETER_SECURITY_POLICY = new JavaParameterSecurityPolicy();

	@Override
	protected void preExecute() throws JobCancelledException {
		super.preExecute();
	}

	@Override
	protected void postExecute() throws JobCancelledException {
		super.postExecute();
	}

	@Override
	protected void cleanUp() {
		super.cleanUp();
	}

	@Override
	protected void cancelRequested() {
		// ignore by default
	}

	public abstract String getSADL();
}
