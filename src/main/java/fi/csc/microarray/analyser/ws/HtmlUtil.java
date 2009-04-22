package fi.csc.microarray.analyser.ws;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import fi.csc.microarray.util.Strings;

public class HtmlUtil {
	
	public static void writeHtmlTable(ResultTableCollector annotations, String[] columns) throws FileNotFoundException {
		String[][] table = annotations.asTable(columns);
		PrintWriter out = new PrintWriter(new FileOutputStream("test.html"));
		out.print("<html>");
		out.print("<table>");
		out.print("<tr>");
		for (int i = 0; i < columns.length; i++) {			
			String name = columns[i].contains(":") ? columns[i].split(":")[1] : columns[i];			
			name = Strings.separateUppercaseChars(Strings.startWithUppercase(name), "-");
			name = name.replaceAll("_", " ");
			out.print("<th>" + name + "</th>");
		}
		out.print("</tr>");
		for (int i = 0; i < table.length; i++) {
			out.print("<tr>");
			for (int j = 0; j < table[i].length; j++) {
				out.print("<td>" + table[i][j] + "</td>");
			}
			out.print("</tr>");
		}
		out.print("</table>");
		out.print("</html>");
		out.close();
	}
	

}
