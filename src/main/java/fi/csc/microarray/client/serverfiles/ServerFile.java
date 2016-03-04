package fi.csc.microarray.client.serverfiles;

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.util.LinkedList;

public class ServerFile extends File {
	
	public static final String SERVER_SESSION_ROOT_FOLDER = "Cloud sessions";

	private String name;
	private String url;
	private LinkedList<ServerFile> children = new LinkedList<ServerFile>();
	private boolean isDirectory;
	private ServerFile parent;
	
	public ServerFile(String url) {
		super(url);
		this.url = url;
		this.isDirectory = url.endsWith("/");
		
		String p = url;
		if (isDirectory) {
			this.name = p.substring(p.substring(0, p.length() - 1).lastIndexOf("/") + 1, p.length() - 1);
		} else {
			this.name = p.substring(p.lastIndexOf("/") + 1, p.length());
		}
	}

	public ServerFile(File dir, String filename) {
		this(concatenate(dir, filename));
	}

	private static String concatenate(File dir, String filename) {
		String dirString = dir.getPath();
		if (dirString.endsWith("/")) {
			dirString = dirString.substring(0, dirString.length() - 1);
		}
		return dirString + "/" + filename;
	}

	public void addChild(ServerFile child) {
		if (!isDirectory) {
			throw new UnsupportedOperationException("cannot add children if not a directory");
		}
		children.add(child);
		child.setParent(this);
	}

	@Override
	public boolean canExecute() {
		return false;
	}

	@Override
	public boolean canWrite() {
		return false;
	}

	@Override
	public int compareTo(File pathname) {
		if (pathname instanceof ServerFile) {
			return url.compareTo(((ServerFile)pathname).url);
		} else {
			throw new IllegalArgumentException();
		}
	}

	@Override
	public boolean canRead() {
		return true;
	};

	@Override
	public boolean createNewFile() throws IOException {
		return false;
	}

	@Override
	public boolean delete() {
		return false;
	}

	@Override
	public void deleteOnExit() {
		// ignore
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ServerFile) {
			return url.equals(((ServerFile)obj).url);
		} else {
			return false;
		}
	};

	@Override
	public boolean exists() {
		return true;
	}

	@Override
	public File getAbsoluteFile() {
		return this;
	}

	@Override
	public String getAbsolutePath() {
		return url;
	}

	@Override
	public File getCanonicalFile() throws IOException {
		return this;
	};

	@Override
	public String getCanonicalPath() throws IOException {
		return url;
	}

	@Override
	public long getFreeSpace() {
		return 0; // cannot write anyway
	}

	@Override
	public String getName() {
		return name;
	};

	@Override
	public String getParent() {
		return super.getParent();
	}

	@Override
	public File getParentFile() {
		return parent;
	}
	
	public void setParent(ServerFile parent) {
		this.parent = parent;
	}

	@Override
	public String getPath() {
		return url;
	}

	@Override
	public long getTotalSpace() {
		return 0; // could calculate this, but is not needed
	}

	@Override
	public long getUsableSpace() {
		return 0; // cannot write anyway
	}

	@Override
	public int hashCode() {
		return url.hashCode();
	};

	@Override
	public boolean isAbsolute() {
		return url.startsWith(SERVER_SESSION_ROOT_FOLDER);
	}

	@Override
	public boolean isDirectory() {
		return isDirectory;
	};

	@Override
	public boolean isFile() {
		return !isDirectory;
	}

	@Override
	public boolean isHidden() {
		return false;
	}

	@Override
	public long lastModified() {
		return 0; // not supported  
	}

	@Override
	public long length() {
		return 0; // not supported
	};

	@Override
	public String[] list() {
		return list(null);
	}

	@Override
	public String[] list(java.io.FilenameFilter filter) {

		LinkedList<String> names = new LinkedList<String>();
		for (ServerFile child : children) {
			if (filter == null || filter.accept(this, child.getName())) {
				names.add(child.getName());
			}
		}
		
		return names.toArray(new String[0]);
	};

	@Override
	public File[] listFiles() {
		return children.toArray(new File[0]);
	};

	@Override
	public File[] listFiles(FileFilter filter) {
		LinkedList<ServerFile> filtered = new LinkedList<ServerFile>();
		
		for (ServerFile child : children) {
			if (filter.accept(child)) {
				filtered.add(child);
			}
		}
		return filtered.toArray(new File[0]);
	}

	@Override
	public File[] listFiles(FilenameFilter filter) {
		LinkedList<ServerFile> filtered = new LinkedList<ServerFile>();
		
		for (ServerFile child : children) {
			if (filter.accept(this, child.getName())) {
				filtered.add(child);
			}
		}
		return filtered.toArray(new File[0]);
	}

	@Override
	public boolean mkdir() {
		return false;
	}

	@Override
	public boolean mkdirs() {
		return false;
	}

	@Override
	public boolean renameTo(File dest) {
		return false;
	}

	@Override
	public boolean setExecutable(boolean executable) {
		return false;
	}

	@Override
	public boolean setExecutable(boolean executable, boolean ownerOnly) {
		return false;
	}

	@Override
	public boolean setLastModified(long time) {
		return false;
	}

	@Override
	public boolean setReadable(boolean readable) {
		return false;
	}

	@Override
	public boolean setReadable(boolean readable, boolean ownerOnly) {
		return false;
	}

	@Override
	public boolean setReadOnly() {
		return false;
	}

	@Override
	public boolean setWritable(boolean writable) {
		return false;
	}

	@Override
	public boolean setWritable(boolean writable, boolean ownerOnly) {
		return false;
	}

	@Override
	public String toString() {
		return getName();
	}

	@Override
	public URI toURI() {
		try {
			return toURL().toURI();
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public URL toURL() throws MalformedURLException {
		return new URL(url);
	}
}
