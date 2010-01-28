package fi.csc.microarray.client.wizard.affymetrix;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;

import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.constants.VisualConstants;


public class FileGroupPanel extends TitlePanel {

    private JTextArea fileList = new JTextArea(7, 30);
	private JComboBox groupList = new JComboBox();
	private SelectionSet selectionSet = new SelectionSet();
	private Map<String, FileGroup> groups = new HashMap<String, FileGroup>();
	private String groupname = null;
	private JTextArea infoArea = new JTextArea(5, 35);
	private FileGroupPanelDesc desc;
	private JFileChooser fileChooser;

	public FileGroupPanel(FileGroupPanelDesc d) {

		super();

		desc = d;
		setTitle("Wizard: Load data files");

		JPanel contentPanel = getContentPanel();

		addContent(contentPanel);
	}

	private JPanel getContentPanel() {
		selectionSet.add("Group1", "filegroup1");
		selectionSet.add("Group2", "filegroup2");
		selectionSet.add("Group3", "filegroup3");
		selectionSet.add("Group4", "filegroup4");
		selectionSet.add("Group5", "filegroup5");
		selectionSet.initComboBox(groupList);
		
		JPanel cPanel = new JPanel();
		cPanel.setLayout(new BorderLayout());
		cPanel.setBorder(new EmptyBorder(new Insets(10, 10, 10, 10)));
		
		JPanel groupPanel = new JPanel();
		groupPanel.setBorder(new EtchedBorder());
		groupPanel.setLayout(new BoxLayout(groupPanel, BoxLayout.Y_AXIS));
		cPanel.add(groupPanel, BorderLayout.NORTH);
		
		JPanel upperPanel = new JPanel();
		upperPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		upperPanel.add(new JLabel("Select group:"));
		upperPanel.add(groupList);
		groupList.addActionListener(new GroupSelector());

		JScrollPane scrollPane = new JScrollPane(fileList,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		fileList.setEditable(false);

		JPanel middlePanel = new JPanel();
		middlePanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		middlePanel.add(scrollPane);

		groupPanel.add(upperPanel);
		groupPanel.add(middlePanel);

		JPanel lowerPanel = new JPanel();
		lowerPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		groupPanel.add(lowerPanel);
		
		JButton addButton = new JButton("Add Files");
		addButton.addActionListener(new FileAdder());

		JButton resetButton = new JButton("Reset group");
		resetButton.addActionListener(new GroupReset());

		lowerPanel.add(addButton);
		lowerPanel.add(resetButton);

		JPanel infoPane = new JPanel();
		infoPane.setLayout(new BoxLayout(infoPane, BoxLayout.Y_AXIS));
		
		JPanel flowPane4 = new JPanel();
		flowPane4.setLayout(new FlowLayout(FlowLayout.LEFT));
		infoPane.add(flowPane4);
		
		infoArea.setEditable(false);
		infoArea.setFont(new Font("Verdana", Font.BOLD, 10));
		//infoArea.setBorder(new EtchedBorder());
		infoArea.setBackground(VisualConstants.TEXTAREA_UNEDITABLE_BACKGROUND);
		
		//JScrollPane infoScrollPane = new JScrollPane(infoArea,
		//		JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
		//		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		//flowPane4.add(infoScrollPane);
		
		flowPane4.add(infoArea);

		cPanel.add(infoPane, BorderLayout.CENTER);

		JPanel textPanel = new JPanel();
		cPanel.add(textPanel, BorderLayout.SOUTH);
		textPanel.setLayout(new BoxLayout(textPanel, BoxLayout.Y_AXIS));
		
		textPanel.add(new JLabel("Specify your experimental design by importing the chips as one or more"));
		textPanel.add(new JLabel("groups. A group can denote a treatment, or some other sample status that is"));
		textPanel.add(new JLabel("of interest. For example, if your experiment consists of healthy controls"));
		textPanel.add(new JLabel("and cancer patients, you should import controls as a group 1, and patients"));
		textPanel.add(new JLabel("as group 2. You can specify a maximum of 5 groups."));
		
		groupname = "filegroup1";

		infoArea.setText(getInfoString());

		return cPanel;
	}

	public Vector<FileGroup> getGroups() {
		FileGroup group;
		Vector<FileGroup> groupSet = new Vector<FileGroup>();
		Iterator<SelectionItem> iter = selectionSet.getIterator();
		
		while (iter.hasNext()) {
			SelectionItem item = iter.next();
			String name = item.getInternalValue();
			
			group = groups.get(name);
			if (group != null) {
				if (group.getNum() > 0)
					groupSet.add(group);
			}
		}
		return groupSet;
	}

	public boolean isGroupsEmpty() {
		FileGroup group;
		Iterator<SelectionItem> iter = selectionSet.getIterator();
		boolean empty = true;

		while (iter.hasNext()) {
			SelectionItem item = iter.next();
			String name = item.getInternalValue();
			
			group = groups.get(name);
			if (group != null)
				if (group.getNum() > 0)
					empty = false;
		}

		return empty;
	}

	private String getInfoString() {
		FileGroup group;
		Iterator<SelectionItem> iter = selectionSet.getIterator();
		String text = null;
		String aux = null;

		while (iter.hasNext()) {
			SelectionItem item = iter.next();
			String publicName = item.getPublicValue();
			String internalName = item.getInternalValue();
			
			group = groups.get(internalName);
			if (group != null)
				if (group.getNum() > 0)
					aux = publicName + ": " + group.getNum() + " files";

			if (aux != null) {
				if (text != null)
					text = text + "\n" + aux;
				else
					text = aux;
			}
			aux = null;
		}

		if (text == null)
			text = new String(" ");

		return text;
	}

	private class GroupSelector implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			String text = null;
			FileGroup group = null;

			Object obj = groupList.getSelectedItem();
			groupname = selectionSet.getInternalValue(obj.toString());

			group = groups.get(groupname);
			if (group != null)
				text = group.getFileNameList();

			fileList.setText(text);
		}
	}

	private class FileAdder implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			FileGroup group = null;
			String text = null;

			if(fileChooser == null){
				fileChooser = ImportUtils.getFixedFileChooser();			
				fileChooser.setMultiSelectionEnabled(true);
			}
			
			int returnVal = fileChooser.showOpenDialog(null);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File[] files = fileChooser.getSelectedFiles();

				group = groups.get(groupname);
				if (group == null) {
					group = new FileGroup(groupname);
					groups.put(groupname, group);
				}

				for (File file : files)
					group.addFile(file);

				text = group.getFileNameList();
				fileList.setText(text);
				infoArea.setText(getInfoString());

				desc.setNextButton();
			}
		}
	}

	private class GroupReset implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			FileGroup group = null;

			group = groups.get(groupname);
			if (group != null)
				group.clear();

			fileList.setText("");
			infoArea.setText(getInfoString());

			desc.setNextButton();
		}
	}
}
