package fi.csc.microarray.client.serverfiles;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.HashMap;
import java.util.List;

import javax.swing.filechooser.FileSystemView;

import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.DerbyMetadataServer;

/**
 * This class provides a view into a generic directory service.
 */
public class ServerFileSystemView extends FileSystemView {

	private ServerFile rootFile;
	
	public ServerFileSystemView(ServerFile rootFile) {
		this.rootFile = rootFile;
	}

	public static ServerFileSystemView parseFromPaths(String prefix, List<DbSession> sessions) throws MalformedURLException {
		
		HashMap<String, ServerFile> dirs = new HashMap<String, ServerFile>();
		ServerFile root = new ServerFile(prefix + "/");
		dirs.put(prefix + "/", root);
		
		for (DbSession session : sessions) {
			String path = session.getName();
			String fullPath = prefix + "/" + path;
			ServerFile file = new ServerFile(fullPath);
			ServerFile parent = dirs.get(fullPath.substring(0, fullPath.substring(0, fullPath.length() - 1).lastIndexOf("/") + 1));
			parent.addChild(file);
			if (file.isDirectory()) {
				dirs.put(fullPath, file);
			}
		}
		
		return new ServerFileSystemView(root);
	}
	
	
	@Override
	public File createNewFolder(File aContainingDir) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public File[] getRoots() {
		return new File[] { rootFile };
	}

	public File getRoot() {
		return rootFile;
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

	public File getExampleSessionDir() {
		for (File file : rootFile.listFiles()) {
			if (file instanceof ServerFile) {
				ServerFile serverFile = (ServerFile) file;
				if (serverFile.getName().equals(DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_FOLDER)) {
					return serverFile;
				}
			}
		}
		return null;
	}
}