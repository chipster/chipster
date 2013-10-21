package fi.csc.chipster.web.tooledit;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringReader;

import com.vaadin.ui.TextArea;
import com.vaadin.ui.VerticalLayout;

import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.analyser.SADLTool.ParsedScript;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.Strings;

/**
 * Text editor
 * @author Gintare Pacauskaite
 *
 */
public class TextEditor extends VerticalLayout{
	private static final long serialVersionUID = -7074541336842177583L;
	
	private static final String HEADER_PREFIX= "#";
	
	private static class SeparatedContent {
		public String header = "";
		public String code = "";
	}

	private TextArea txtArea;

	
	public TextEditor(ToolEditorUI root) {
		init();
	}
	
	private void init() {
		txtArea = new TextArea();
		txtArea.setSizeFull();
		// for some reasons size full does not do anything to height
		txtArea.setRows(50);
		
		this.addComponent(txtArea);
	}
	
	
	/**
	 * 
	 * @return plain header without comment prefix
	 * @throws IOException 
	 */
	public String getHeader() throws IOException {
		SADLTool tool = new SADLTool("#");

		// this will make sure that there's empty line after header and non empty line in the beginning of the code
		SeparatedContent content = parseContent();
		
		// parse header to sadl
		ParsedScript parsedScript = tool.parseScript(new ByteArrayInputStream(content.header.getBytes()));
		
		// update text area (now with exactly one empty line between header and code
		txtArea.setValue(content.header + content.code);
		
		return parsedScript.SADL;
	}

	
	/**
	 * 
	 * @param newHeader without prefixes
	 * @throws IOException
	 */
	public void setHeader(String newHeader) throws IOException {
		
		// save old content
		// parsed code will start with non empty line
		SeparatedContent oldContent = parseContent();
		
		// prefix new header, add empty line in the end
		SADLTool sadlTool = new SADLTool("#");
		ParsedScript temp = new ParsedScript();
		temp.SADL = newHeader;
		String prefixedHeader = sadlTool.toScriptString(temp) + "\n";
		
		// set new content
		txtArea.setValue(prefixedHeader + oldContent.code);
		
	}
	
	
	public void clearAllText() {
		txtArea.setValue("");
	}

	/**
	 * Resulting header will end in one empty line and resulting code will start with non empty line.
	 * @return
	 * @throws IOException
	 */
	private SeparatedContent parseContent() throws IOException {
		BufferedReader reader = new BufferedReader(new StringReader(txtArea.getValue()));
		StringBuilder header = new StringBuilder();
		StringBuilder code = new StringBuilder();
		boolean parsingHeader = true;
		boolean firstSignificantLineRead = false;
		try {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {

				// eat empty lines and empty prefixed lines from beginning
				if (!firstSignificantLineRead && (
						(line.trim().isEmpty()) ||
						(line.startsWith(HEADER_PREFIX) && line.substring(HEADER_PREFIX.length()).trim().isEmpty())
						)) {
					continue;
				}

				// first significant line
				if (!firstSignificantLineRead) {
					// first line is prefixed -> header
					if (line.startsWith(HEADER_PREFIX)) {
						parsingHeader = true;
						header.append(line + "\n");
					} 

					// first line not prefixed -> not header
					else {
						parsingHeader = false;
						code.append(line + "\n");
					}

					firstSignificantLineRead = true;
					continue;
				}

				// rest of the lines, header ends with empty or non-prefixed line or prefixed white space line
				if (parsingHeader && line.startsWith(HEADER_PREFIX) && !line.substring(HEADER_PREFIX.length()).trim().isEmpty()) {
					header.append(line + "\n");
				} else {
					code.append(line + "\n");
					parsingHeader = false;
				}
			}
		} finally {
			IOUtils.closeIfPossible(reader);
		}
		
		// add empty line to the end of header, remove empty lines from the beginning of code
		SeparatedContent content = new SeparatedContent();
		content.header = header.toString() + "\n";
		content.code = Strings.removeEmptyLinesFromBeginning(code.toString());
		return content;
	}

}
