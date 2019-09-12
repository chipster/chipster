package fi.csc.microarray.messaging.message;

import java.util.Iterator;
import java.util.List;

import fi.csc.microarray.comp.ToolDescription;
import fi.csc.microarray.comp.ToolDescription.ParameterDescription;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;

public class JobMessageUtils {
	/**
	 * This should really be in the GenericJobMessage, but static methods in
	 * interfaces are only allowed starting from Java 1.8.
	 * 
	 * @param securityPolicy
	 * @param description
	 * @param parameters
	 * @return
	 * @throws ParameterValidityException
	 */
	public static List<String> checkParameterSafety(ParameterSecurityPolicy securityPolicy, ToolDescription description,
			List<String> parameters) throws ParameterValidityException {
		// Do argument checking first
		if (securityPolicy == null) {
			throw new IllegalArgumentException("security policy cannot be null");
		}
		if (description == null) {
			throw new IllegalArgumentException("tool description cannot be null");
		}

		// Count parameter descriptions
		int parameterDescriptionCount = 0;
		for (Iterator<ParameterDescription> iterator = description.getParameters().iterator(); iterator
				.hasNext(); iterator.next()) {
			parameterDescriptionCount++;
		}

		// Check that description and values match
		if (parameterDescriptionCount != parameters.size()) {
			throw new IllegalArgumentException(
					"number of parameter descriptions does not match the number of parameter values");
		}

		// Validate parameters
		Iterator<ParameterDescription> descriptionIterator = description.getParameters().iterator();
		for (String parameter : parameters) {
			ParameterDescription parameterDescription = descriptionIterator.next();

			if (parameterDescription.isChecked()) {
				if (!securityPolicy.isValueValid(parameter, parameterDescription)) {
					throw new ParameterValidityException(
							"illegal value for parameter " + parameterDescription.getName() + ": " + parameter);
				}
			} else {
				if (!securityPolicy.allowUncheckedParameters(description)) {
					throw new UnsupportedOperationException("unchecked parameters are not allowed");
				}
			}
		}

		// Everything was ok, return the parameters
		return parameters;
	}
}
