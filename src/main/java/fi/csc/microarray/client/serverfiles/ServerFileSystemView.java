package fi.csc.microarray.client.serverfiles;

import java.io.File;
import java.io.IOException;

import javax.swing.filechooser.FileSystemView;

/**
 * This class provides a view into a generic directory service.
 */
public class ServerFileSystemView extends FileSystemView {

	private ServerFile rootFile;
	
	public ServerFileSystemView(ServerFile rootFile) {
		this.rootFile = rootFile;
	}
	
	@Override
	public File createNewFolder(File aContainingDir) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public File[] getRoots() {
		return new File[] { rootFile };
	}

	@Override
	public boolean isHiddenFile(File f) {
		return false;
	}

	@Override
	public boolean isRoot(File f) {
		return rootFile.equals(f);
	}

	@Override
	public File getHomeDirectory() {
		return rootFile;
	}

	@Override
	public File[] getFiles(File dir, boolean useFileHiding) {
		return dir.listFiles();
	}

	@Override
	public String getSystemDisplayName(File f) {
		return f.getName();
	};
	
	@Override
	public File getParentDirectory(File dir) {
		return dir.getParentFile();
	}
	
	@Override
	public boolean isDrive(File dir) {
		return rootFile.equals(dir);
	}
}