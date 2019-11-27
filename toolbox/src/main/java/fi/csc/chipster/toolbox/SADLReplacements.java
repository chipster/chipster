package fi.csc.chipster.toolbox;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import fi.csc.microarray.description.SADLDescription.Name;

public class SADLReplacements {

	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(SADLReplacements.class);

	public static final String FILES = "FILES";
	public static final String SYMLINK_TARGET = "SYMLINK_TARGET";

	private File toolsBin;

	public SADLReplacements(File toolsBin) {
		this.toolsBin = toolsBin;
	}

	public Name[] processNames(Collection<Name> options) throws IOException {
		ArrayList<Name> newOptions = new ArrayList<>();
		for (Name option : options) {
			List<String> replaced = processReplacements(option.getID());
			if (replaced.size() == 1 && option.getID().equals(replaced.get(0))) {
				// there is no replacement, keep the original Name object
				newOptions.add(option);
			} else {
				// create Name objects for all replacement strings
				newOptions.addAll(getNames(replaced));
			}
		}
		return newOptions.toArray(new Name[0]);
	}

	public String[] processStrings(Collection<String> options) throws IOException {
		ArrayList<String> newOptions = new ArrayList<>();
		for (String option : options) {
			return processReplacements(option).toArray(new String[0]);
		}
		return newOptions.toArray(new String[0]);
	}

	private List<String> processReplacements(String option) throws IOException {

		// split by white space
		String[] tokens = parseReplacements(option);
		if (tokens.length >= 1) {
			if (FILES.equals(tokens[0])) {
				if (tokens.length == 3) {
					return getFiles(tokens[1], tokens[2]);
				} else {
					throw new IOException(FILES + " expected 2 parameters, but got " + (tokens.length - 1));
				}
			} else if (SYMLINK_TARGET.equals(tokens[0])) {
				if (tokens.length == 3) {
					return Arrays.asList(new String[] { getSymlinkTarget(tokens[1], tokens[2]) });
				} else {
					throw new IOException(SYMLINK_TARGET + " expected 2 parameters, but got " + (tokens.length - 1));
				}
			}
		}
		return Arrays.asList(new String[] { option });
	}

	private String[] parseReplacements(String value) {
		// split by white space
		String[] tokens = value.split("\\s+");
		return tokens;
	}

	private ArrayList<Name> getNames(Collection<String> ids) {
		ArrayList<Name> names = new ArrayList<>();
		for (String id : ids) {
			Name name = new Name();
			name.setID(id);
			names.add(name);
		}
		return names;
	}

	private ArrayList<String> getFiles(String path, String endsWith) throws IOException {
		ArrayList<String> basenames = new ArrayList<>();
		File dir = new File(toolsBin, path);
		if (!dir.exists()) {
			throw new IOException("failed to get enum options from files, file not found: " + dir.getPath());
		}
		for (String fileName : dir.list()) {
			if (fileName.endsWith(endsWith)) {
				String basename = fileName.substring(0, fileName.length() - endsWith.length());
				basenames.add(basename);
			}
		}

		// sort alphabetically
		basenames.sort(null);

		return basenames;
	}

	private String getSymlinkTarget(String path, String fileExtension) throws IOException {
		Path symlink = new File(toolsBin, path).toPath();
		if (Files.isSymbolicLink(symlink)) {
			String target = Files.readSymbolicLink(symlink).getFileName().toString();
			target = StringUtils.removeEnd(target, fileExtension);
			return target;
		} else {
			throw new IOException("failed to get symlink target, not a symlink: " + path);
		}
	}
}
