package fi.csc.microarray.analyser.ws;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.Strings;

public class HtmlUtil {
	
	public static interface ValueHtmlFormatter {
		public String format(String value, String[] currentRow);		
	}
	
	public static ValueHtmlFormatter NO_FORMATTING_FORMATTER = new ValueHtmlFormatter() {
		public String format(String value, String[] currentRow) {
			return value;
		}		
	};
	
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
			out.println("<html>");
			out.println("<head>");
			out.println(css);
			out.println("<title>" + title + "</title>");
			out.println("</head>");			
		}
		
		public void closeDocument() {
			out.println("</html>");
		}
		
		public void close() {
			out.close();
		}

		public void openTable() {
			out.println("<table>");			
		}

		public void closeTable() {
			out.println("</table>");			
		}

		public void openRow() {
			out.println("<tr>");
		}
		
		public void closeRow() {
			out.println("</tr>");
		}

		public void headerCell(String name) {
			out.println("<th>" + name + "</th>");

		}

		public void cell(String data) {
			out.println("<td>" + data + "</td>");
		}

		public void heading(String heading) {
			out.println("<h1>" + heading + "</h1>");			
		}
	}
	
	public static void writeHtmlTable(ResultTableCollector annotations, String[] columnIds, String[] columnNames, ValueHtmlFormatter[] valueHtmlFormatters, String title, File output) throws FileNotFoundException {
		FileOutputStream out = null;
		try {
			out = new FileOutputStream(output);
			writeHtmlTable(annotations, columnIds, columnNames, valueHtmlFormatters, title, out);
			
		} finally {
			IOUtils.closeIfPossible(out);
		}
		
	}
	public static void writeHtmlTable(ResultTableCollector annotations, String[] columnIds, String[] columnNames, ValueHtmlFormatter[] valueHtmlFormatters, String title, OutputStream outputStream) throws FileNotFoundException {
		HtmlTemplate out = new HtmlTemplate(new PrintWriter(outputStream));
		out.openDocument(title);
		
		// write header row
		out.heading(title);
		out.openTable();
		out.openRow();
		for (String columnName : columnNames) {			
			out.headerCell(columnName);
		}
		out.closeRow();

		// write data rows
		String[][] table = annotations.asTable(columnIds);
		for (int i = 0; i < table.length; i++) {
			out.openRow();
			for (int j = 0; j < table[i].length; j++) {
				String string = table[i][j];
				if (valueHtmlFormatters != null) {
					string = valueHtmlFormatters[j].format(string, table[i]);
				}
				out.cell(string);
			}
			out.closeRow();
		}
		out.closeTable();
		out.closeDocument();
		out.close();
	}
	
	public static void writeTextTable(ResultTableCollector annotations, String[] columnIds, String[] columnNames, OutputStream outputStream) throws FileNotFoundException {
		PrintWriter out = new PrintWriter(new OutputStreamWriter(outputStream));
		String delimeter = "\t";

		// write header row
		out.print(Strings.delimit(Arrays.asList(columnNames), delimeter));		
		out.print("\n");

		// write data rows
		String[][] table = annotations.asTable(columnIds);
		for (int i = 0; i < table.length; i++) {
			out.print(Strings.delimit(Arrays.asList(table[i]), delimeter));
			out.print("\n");
		}
		out.close();
	}

}
