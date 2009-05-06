package fi.csc.microarray.analyser.ws;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import fi.csc.microarray.util.Strings;

public class HtmlUtil {
	
	private static class HtmlTemplate {

		final String css = "<style type=\"text/css\">\n" + 
		"table {\n" + 
		"	border-width: 1px 1px 1px 1px;\n" + 
		"	border-spacing: 2px;\n" + 
		"	border-style: outset outset outset outset;\n" + 
		"	border-color: black black black black;\n" + 
		"	border-collapse: separate;\n" + 
		"	background-color: white;\n" + 
		"}\n" + 
		"table th {\n" + 
		"	border-width: 1px 1px 1px 1px;\n" + 
		"	padding: 1px 1px 1px 1px;\n" + 
		"	border-style: inset inset inset inset;\n" + 
		"	border-color: gray gray gray gray;\n" + 
		"	background-color: white;\n" + 
		"	-moz-border-radius: 0px 0px 0px 0px;\n" + 
		"}\n" + 
		"table td {\n" + 
		"	border-width: 1px 1px 1px 1px;\n" + 
		"	padding: 1px 1px 1px 1px;\n" + 
		"	border-style: inset inset inset inset;\n" + 
		"	border-color: gray gray gray gray;\n" + 
		"	background-color: white;\n" + 
		"	-moz-border-radius: 0px 0px 0px 0px;\n" + 
		"}\n" + 
		"</style>\n";

		PrintWriter out;

		public HtmlTemplate(PrintWriter out) {
			this.out = out;
		}

		public void openDocument(String title) {
			out.print("<html>");
			out.print("<head>");
			out.print(css);
			out.print("<title>" + title + "</title>");
			out.print("</head>");			
		}
		
		public void closeDocument() {
			out.print("</html>");
		}
		
		public void close() {
			out.close();
		}

		public void openTable() {
			out.print("<table>");			
		}

		public void closeTable() {
			out.print("</table>");			
		}

		public void openRow() {
			out.print("<tr>");
		}
		
		public void closeRow() {
			out.print("</tr>");
		}

		public void headerCell(String name) {
			out.print("<th>" + name + "</th>");

		}

		public void cell(String data) {
			out.print("<td>" + data + "</td>");
		}

		public void heading(String heading) {
			out.print("<h1>" + heading + "</h1>");			
		}
	}
	
	public static void writeHtmlTable(ResultTableCollector annotations, String[] columns, String title, File file) throws FileNotFoundException {
		String[][] table = annotations.asTable(columns);
		HtmlTemplate out = new HtmlTemplate(new PrintWriter(new FileOutputStream(file)));
		out.openDocument(title);
		out.heading(title);
		out.openTable();
		out.openRow();
		for (int i = 0; i < columns.length; i++) {			
			String name = columns[i].contains(":") ? columns[i].split(":")[1] : columns[i];			
			name = Strings.separateUppercaseChars(Strings.startWithUppercase(name), "-");
			name = name.replaceAll("_", " ");
			out.headerCell(name);
		}
		out.closeRow();
		for (int i = 0; i < table.length; i++) {
			out.openRow();
			for (int j = 0; j < table[i].length; j++) {
				out.cell(table[i][j]);
			}
			out.closeRow();
		}
		out.closeTable();
		out.closeDocument();
		out.close();
	}
	

}
