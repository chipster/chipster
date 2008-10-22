package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

public class AuthenticationMessage extends CommandMessage {
	public enum AuthenticationOperation {
		LOGIN,
		LOGOUT,
		REQUEST,
		LOGIN_FAILED,
		LOGIN_SUCCEEDED;
		
		public String toString() {
			switch (this) {
			case LOGIN:
				return "login";
			case LOGOUT:
				return "logout";
			case REQUEST:
				return "request";
			case LOGIN_FAILED:
				return "login-failed";
			case LOGIN_SUCCEEDED:
				return "login-succeeded";
			default:
				throw new IllegalArgumentException("unknown operation: " + this);	
			}
		}
	}

	public static final String KEY_PASSWORD = "password";
	
	private String password;

	public AuthenticationMessage() {
		super();
	}
	
	public AuthenticationMessage(AuthenticationOperation operation) {
		super(operation.toString());
	}

	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.password = from.getString(KEY_PASSWORD);
	}

	@Override
	public void marshal(MapMessage to) throws JMSException {
		super.marshal(to);
		to.setString(KEY_PASSWORD, password);
	}
	
	public boolean isRequestForAuthentication() {
		return AuthenticationOperation.REQUEST.toString().equals(getCommand()) || 
			AuthenticationOperation.LOGIN_FAILED.toString().equals(getCommand());
	}
	
	public boolean isLogin() {
		return AuthenticationOperation.LOGIN.toString().equals(getCommand());
	}
	
	public boolean isLogout() {
		return AuthenticationOperation.LOGOUT.toString().equals(getCommand());
	}

	public boolean isLoginAck() {
		return AuthenticationOperation.LOGIN_SUCCEEDED.toString().equals(getCommand());
	}
	
	/**
	 * Returns the password part of credentials associated with this message.
	 */
	public String getPassword() {
		return password;
	}

	/**
	 * @see #getPassword()
	 */
	public void setPassword(String password) {
		this.password = password;
	}

}
