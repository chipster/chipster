package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.Table;

@Entity
@Table(name="accounts")
public class AccountEntry implements Serializable {

	@Id
	@Column(name="username")
	private String username;
	@Column(name="ignoreInStatistics")
	private boolean ignoreInStatiscs;
	
	public AccountEntry() {};
	
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = username;
	}
	public boolean getIgnoreInStatiscs() {
		return ignoreInStatiscs;
	}
	public void setIgnoreInStatiscs(boolean ignoreInStatiscs) {
		this.ignoreInStatiscs = ignoreInStatiscs;
	}
}
