package fi.csc.microarray.client.serverfiles;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import javax.swing.filechooser.FileSystemView;

/**
 * This class provides a view into a generic directory service.
 */
public class ServerFileSystemView extends FileSystemView {

	private ServerFile rootFile;
	
	public ServerFileSystemView(ServerFile rootFile) {
		this.rootFile = rootFile;
	}

	public static ServerFileSystemView parseFromPaths(URL[] urls) throws MalformedURLException {
		
		HashMap<String, ServerFile> dirs = new HashMap<String, ServerFile>();
		ServerFile root = null;
		
		for (URL url : urls) {
			ServerFile file = new ServerFile(url);
			String p = url.getPath();
			ServerFile parent = dirs.get(p.substring(0, p.substring(0, p.length() - 1).lastIndexOf("/") + 1));
			if (parent == null) {
				root = file; // this was root
			} else {
				parent.addChild(file);
			}
			if (file.isDirectory()) {
				dirs.put(url.getPath(), file);
			}
		}
//		
//		ServerFile root = new ServerFile(new URL(prefix + paths[0]));
//		ServerFile home = new ServerFile(new URL(prefix + paths[1]));
//		root.addChild(home);
//		ServerFile bam1 = new ServerFile(new URL(prefix + paths[2]));
//		ServerFile bam2 = new ServerFile(new URL(prefix + paths[3]));
//		home.addChild(bam1);
//		home.addChild(bam2);
		
		if (root == null) {
			throw new IllegalArgumentException("paths were missing root");
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
}