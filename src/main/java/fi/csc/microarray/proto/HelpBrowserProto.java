package fi.csc.microarray.proto;

import java.io.File;

import javax.swing.JFrame;

import org.xhtmlrenderer.simple.FSScrollPane;
import org.xhtmlrenderer.simple.XHTMLPanel;

public class HelpBrowserProto {

	public static void main(String[] args) throws Exception {

		XHTMLPanel panel = new XHTMLPanel();
		panel.setDocument(new File("data/help/microarray-manual/index_fixed_xhtml.html"));
		
		FSScrollPane scroll = new FSScrollPane(panel);
		JFrame frame = new JFrame("test");
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		frame.getContentPane().add(scroll);
		frame.pack();
		frame.setSize(1024, 768);
		frame.setVisible(true);
		
		
	}
}
