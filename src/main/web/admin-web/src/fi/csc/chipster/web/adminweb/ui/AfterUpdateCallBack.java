package fi.csc.chipster.web.adminweb.ui;

public interface AfterUpdateCallBack {
	/**
	 * Hook for things to do after the update. Safe for UI changes.
	 */
	public void updateDone();
}
