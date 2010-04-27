package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.ButtonModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Row;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;


/**
 * @author Petri Klemelä
 */
public class GenomeBrowser extends Visualisation implements ActionListener {

	private static final int CHROMOSOME_COUNT = 22;

	public GenomeBrowser(VisualisationFrame frame) {
		super(frame);
	}

	private static final String ANNOTATION_PATH = "http://chipster-devel.csc.fi:8080/public/ensembl/";
	private static final String CONTENTS_FILE = "contents.txt";
	private static final boolean DEBUG_MODE = true;

	protected JPanel paramPanel;
	private JPanel settingsPanel;

	private JButton useButton;

	protected final ClientApplication application = Session.getSession().getApplication();

	protected GenomePlot plot;

	protected DataBean data;
	private ButtonGroup views;
	private JTextArea megaLocation;
	private JTextArea kiloLocation;
	private JTextArea unitLocation;
	private JComboBox chrBox;
	private JRadioButton horizView;
	private JRadioButton circularView;
	private List<JCheckBox> trackBoxes = new ArrayList<JCheckBox>();

	@Override
	public JPanel getParameterPanel() {

		if (paramPanel == null || data != application.getSelectionManager().getSelectedDataBean()) {

			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

			JPanel settings = this.createSettingsPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			c.weightx = 1.0;
			c.insets.set(5, 0, 0, 0);

			paramPanel.add(tabPane, c);

		}

		return paramPanel;
	}

	public void createAnnotationComponents(JPanel panel, GridBagConstraints c) {

		try {

			InputStream contentsStream;
			if (DEBUG_MODE) {
				contentsStream = new FileInputStream(new File(GenomePlot.FILE_ROOT, "annotations/contents.txt"));
			} else {
				contentsStream = new URL(ANNOTATION_PATH + CONTENTS_FILE).openStream();
			}

			
			List<Row> contents = new AnnotationContents().getContents(contentsStream);

			boolean first = true;
			for (Row row : contents) {

				// FIXME Check if there are other species or versions
				if (first) {
					first = false;
					c.gridy++;
					panel.add(new JLabel(row.species), c);
					c.gridy++;
					panel.add(new JLabel(row.version), c);
				}

				c.gridy++;
				JCheckBox box = new JCheckBox(row.content);
				panel.add(box, c);
				trackBoxes.add(box);
			}

		} catch (IOException e) {
			application.reportException(e);
		}
	}

	public JPanel createSettingsPanel() {

		settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		useButton = new JButton("Draw");
		useButton.addActionListener(this);

		chrBox = new JComboBox();

		// FIXME These should be read from user data file
		for (int i = 1; i <= CHROMOSOME_COUNT; i++) {
			chrBox.addItem(i);
		}

		chrBox.setEnabled(false);

		megaLocation = new JTextArea(1, 3);
		kiloLocation = new JTextArea(1, 3);
		unitLocation = new JTextArea(1, 3);

		JTextArea megaLabel = new JTextArea("M");
		megaLabel.setEditable(false);
		JTextArea kiloLabel = new JTextArea("k");
		kiloLabel.setEditable(false);

		horizView = new JRadioButton("Horizontal");
		horizView.setSelected(true);
		circularView = new JRadioButton("Circular");

		views = new ButtonGroup();
		views.add(horizView);
		views.add(circularView);

		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.gridx = 0;
		c.insets.set(5, 10, 5, 10);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;
		c.gridx = 0;
		c.gridwidth = 5;

		settingsPanel.add(new JLabel("Chromosome"), c);
		c.gridy++;
		settingsPanel.add(chrBox, c);
		c.gridy++;
		// settingsPanel.add(new JLabel("Location"),c);
		//		
		// c.gridy++;
		// c.gridwidth = 1;
		// c.insets.set(5, 10, 5, 0);
		// settingsPanel.add(megaLocation, c);
		// c.gridx++;
		// c.insets.set(5, 0, 5, 0);
		// settingsPanel.add(megaLabel, c);
		// c.gridx++;
		// settingsPanel.add(kiloLocation, c);
		// c.gridx++;
		// settingsPanel.add(kiloLabel, c);
		// c.gridx++;
		// c.insets.set(5, 0, 5, 10);
		// settingsPanel.add(unitLocation, c);

		c.gridx = 0;
		c.gridwidth = 5;
		c.insets.set(5, 10, 5, 10);
		createAnnotationComponents(settingsPanel, c);

		c.gridy++;
		settingsPanel.add(new JLabel("View mode"), c);
		c.gridy++;
		settingsPanel.add(horizView, c);
		c.gridy++;
		settingsPanel.add(circularView, c);

		c.gridy++;
		settingsPanel.add(useButton, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		settingsPanel.add(new JPanel(), c);

		return settingsPanel;
	}

	protected JComponent getColorLabel() {
		return new JLabel("Color: ");
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == useButton) {
			useButtonPressed();
		}
	}

	protected void useButtonPressed() {

		List<Variable> vars = new ArrayList<Variable>();
		vars.add(new ObjVariable(chrBox.getSelectedItem()));
		vars.add(new ObjVariable(megaLocation.toString()));
		vars.add(new ObjVariable(kiloLocation.toString()));
		vars.add(new ObjVariable(unitLocation.toString()));
		vars.add(new ObjVariable(views.getSelection()));

		for (JCheckBox box : trackBoxes) {
			vars.add(new ObjVariable(box));
		}

		application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, VisualisationMethod.GBROWSER, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));

	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		GenomePlot plot;

		List<Variable> vars = getFrame().getVariables();
		if (vars == null || vars.size() < 5) {
			plot = new GenomePlot(1, true, false, false, false, false);
			
		} else {

			boolean horizontal = ((ButtonModel) ((ObjVariable) vars.get(4)).obj).equals(horizView.getModel());

			plot = new GenomePlot(1, // (Integer)((ObjVariable)vars.get(0)).obj,
			horizontal, ((JCheckBox) ((ObjVariable) vars.get(5)).obj).isSelected(), ((JCheckBox) ((ObjVariable) vars.get(6)).obj).isSelected(), ((JCheckBox) ((ObjVariable) vars.get(7)).obj).isSelected(), ((JCheckBox) ((ObjVariable) vars.get(8)).obj).isSelected());

		}

		ChartPanel panel = new ChartPanel(new JFreeChart(plot));
		// panel.setPreferredSize(new Dimension(800, 600));
		plot.chartPanel = panel;
		// SelectableChartPanel selPanel = new SelectableChartPanel(new JFreeChart(plot), plot);
		// selPanel.getChartPanel().addChartMouseListener(plot);

		// selPanel.getChartPanel().addMouseWheelListener(plot);

		panel.setCursor(new Cursor(Cursor.HAND_CURSOR));

		// panel.addMouseWheelListener(plot);

		for (View view : plot.getViews()) {
			panel.addMouseListener(view);
			panel.addMouseMotionListener(view);
			panel.addMouseWheelListener(view);
		}

		return panel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return true;
	}
	
	public class ObjVariable extends Variable {

		public Object obj;

		public ObjVariable(Object obj) {
			super(null, null);

			this.obj = obj;
		}
	}
}
