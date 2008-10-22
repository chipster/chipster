package fi.csc.microarray.wizard.affymetrix;

import java.util.*;
import java.io.*;


public class FileGroup {
	private ArrayList<File> list = null;

	private String name;

	private int num = 0;

	public FileGroup(String name) {
		this.name = name;
	}

	public int getNum() {
		return num;
	}

	public String getName() {
		return name;
	}

	public ArrayList<File> getFiles() {
		return list;
	}

	public void addFile(File file) {
		ListIterator<File> iter = null;
		Boolean namefound = false;
		File listFile = null;
		String fileName = file.getName();
		String listName = null;

		if (list == null) {
			list = new ArrayList<File>();

			list.add(file);
			++num;

			return;
		}

		iter = list.listIterator();
		while (iter.hasNext()) {
			listFile = (File) iter.next();
			listName = listFile.getName();
			if (listName.equals(fileName))
				namefound = true;
		}

		if (namefound == false) {
			list.add(file);
			++num;
		}
	}

	public String getFileNameList() {
		String text = null;
		ListIterator<File> iter = null;

		if (list == null)
			return null;

		iter = list.listIterator();
		while (iter.hasNext()) {
			File file = (File) iter.next();
			String name = file.getName();
			//	    String name = file.getAbsolutePath();
			if (text != null) {
				text = text + "\n" + name;
			} else {
				text = name;
			}
		}

		return text;
	}

	public void clear() {
		list.clear();
		list = null;
		num = 0;
	}
}
