package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

import javax.persistence.Column;

public class JobLogEntryRowCount extends JobLogEntry implements Serializable {

	@Column(name=StatDataSource.ROW_COUNT)
	private int rowCount;
	
	public int getRowCount() {
		return rowCount;
	}
	public void setRowCount(int rowCount) {
		this.rowCount = rowCount;
	}
}
