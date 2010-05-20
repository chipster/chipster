package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Row;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.handlers.LocalFileDataBeanHandler;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.util.IOUtils;


/**
 * @author Petri Klemelä, Aleksi Kallio
 */
public class GenomeBrowser extends Visualisation implements ActionListener, RegionListener, VetoableChangeListener {

	private static final String ANNOTATION_URL_PATH = "annotations";

	private static final int CHROMOSOME_COUNT = 22;

	private static final String CONTENTS_FILE = "contents.txt";

	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";

	private JPanel paramPanel;
	private JPanel settingsPanel;
	private JPanel plotPanel = new JPanel(new CardLayout());
		
	private JButton drawButton;

	private final ClientApplication application = Session.getSession().getApplication();

	private GenomePlot plot;

	private DataBean data;
	private JTextArea megaLocation;
	private JTextArea kiloLocation;
	private JTextArea unitLocation;
	private JComboBox chrBox;
	private JComboBox genomeBox;
//	private JRadioButton horizView;
//	private JRadioButton circularView;
	private List<JCheckBox> trackBoxes = new ArrayList<JCheckBox>();

	public GenomeBrowser(VisualisationFrame frame) {
		super(frame);
	}

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

		InputStream contentsStream = null;

		try {
			// parse what annotations we have available
			contentsStream = new URL(fetchAnnotationUrl() + "/" + CONTENTS_FILE).openStream();
			AnnotationContents annotationContentFile = new AnnotationContents();
			annotationContentFile.parseFrom(contentsStream);
			List<Row> contents = annotationContentFile.getRows(); 

			// read genome name and version for each annotation file
			LinkedHashSet<String> genomes = annotationContentFile.getGenomes();
			c.gridy++;
			settingsPanel.add(new JLabel("Genome"), c);			
			c.gridy++;
			genomeBox = new JComboBox();
			for (String genome : genomes) {
				genomeBox.addItem(genome);
			}
			panel.add(genomeBox, c);

			// list available chromosomes
			chrBox = new JComboBox();

			// FIXME These should be read from user data file
			for (int i = 1; i <= CHROMOSOME_COUNT; i++) {
				chrBox.addItem(""+i);
			}

			c.gridy++;
			settingsPanel.add(new JLabel("Chromosome"), c);
			c.gridy++;
			settingsPanel.add(chrBox, c);
			
			// list available track types for the genome
			for (Row row : contents) {
				if (genomeBox.getSelectedItem().equals(row.version)) {
					c.gridy++;
					JCheckBox box = new JCheckBox(row.content);
					panel.add(box, c);
					trackBoxes.add(box);
				}
			}

		} catch (IOException e) {
			application.reportException(e);
			
		} finally {
			IOUtils.closeIfPossible(contentsStream);
		}
	}

	public JPanel createSettingsPanel() {

		settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		drawButton = new JButton("Draw");
		drawButton.addActionListener(this);

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
		
		megaLocation = new JTextArea(1, 3);
		megaLocation.addVetoableChangeListener(this);
		kiloLocation = new JTextArea(1, 3);
		unitLocation = new JTextArea(1, 3);

		JTextArea megaLabel = new JTextArea("M");
		megaLabel.setEditable(false);
		JTextArea kiloLabel = new JTextArea("k");
		kiloLabel.setEditable(false);

		settingsPanel.add(new JLabel("Location"),c);

		c.gridy++;
		c.gridwidth = 1;
		c.insets.set(5, 10, 5, 0);
		settingsPanel.add(megaLocation, c);
		c.gridx++;
		c.insets.set(5, 0, 5, 0);
		settingsPanel.add(megaLabel, c);
		c.gridx++;
		settingsPanel.add(kiloLocation, c);
		c.gridx++;
		settingsPanel.add(kiloLabel, c);
		c.gridx++;
		c.insets.set(5, 0, 5, 10);
		settingsPanel.add(unitLocation, c);

		c.gridx = 0;
		c.gridwidth = 5;
		c.insets.set(5, 10, 5, 10);
		createAnnotationComponents(settingsPanel, c);

//		horizView = new JRadioButton("Horizontal");
//		horizView.setSelected(true);
//		circularView = new JRadioButton("Circular");
//
//		views = new ButtonGroup();
//		views.add(horizView);
//		views.add(circularView);
//
//		c.gridy++;
//		settingsPanel.add(new JLabel("View mode"), c);
//		c.gridy++;
//		settingsPanel.add(horizView, c);
//		c.gridy++;
//		settingsPanel.add(circularView, c);

		c.gridy++;
		settingsPanel.add(drawButton, c);
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

		if (source == drawButton) {
			showVisualisation();
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		this.data = data;
		
		// create panel with card layout and put message panel there
		JPanel waitPanel = new JPanel();
		waitPanel.add(new JLabel("Please select parameters"));
		plotPanel.add(waitPanel, WAITPANEL);
		
		return plotPanel;	
	}

	private void showVisualisation() {

		try {
			// get local data
			LocalFileDataBeanHandler handler = (LocalFileDataBeanHandler)data.getHandler();

			// fetch emote annotation data
			URL annotationUrl = fetchAnnotationUrl();
			String genome = (String)genomeBox.getSelectedItem();

			// create the plot
			GenomePlot plot = new GenomePlot(true);
			TrackFactory.addCytobandTracks(plot, new DataSource(annotationUrl, "Homo_sapiens.GRCh37.57_karyotype.tsv")); // using always the same
			TrackFactory.addGeneTracks(plot, new DataSource(annotationUrl, "Homo_sapiens." + genome + "_genes.tsv"));
			TrackFactory.addReadTracks(plot, new DataSource[] { new DataSource(handler.getFile(data))}, new DataSource(annotationUrl, "Homo_sapiens." + genome + "_seq.tsv"));
			TrackFactory.addRulerTrack(plot);
			plot.start("1", 1024 * 1024 * 250d);
			plot.addDataRegionListener(this);
			
			// wrap it in a panel
			ChartPanel panel = new ChartPanel(new JFreeChart(plot));
			plot.chartPanel = panel;
			panel.setCursor(new Cursor(Cursor.HAND_CURSOR));

			// add mouse listeners
			for (View view : plot.getViews()) {
				panel.addMouseListener(view);
				panel.addMouseMotionListener(view);
				panel.addMouseWheelListener(view);
			}

			// put panel on top of card layout
			if (plotPanel.getComponentCount() == 2) {
				plotPanel.remove(1);
			}
			plotPanel.add(panel, PLOTPANEL);
		    CardLayout cl = (CardLayout)(plotPanel.getLayout());
		    cl.show(plotPanel, PLOTPANEL);

		} catch (Exception e) {
			application.reportException(e);
		}
	}

	private URL fetchAnnotationUrl() {
		try {
			MessagingEndpoint messagingEndpoint = Session.getSession().getMessagingEndpoint("client-endpoint");
			FileBrokerClient fileBrokerClient = new FileBrokerClient(messagingEndpoint.createTopic(Topics.Name.URL_TOPIC, AccessMode.WRITE));
			URL annotationUrl = new URL(fileBrokerClient.getPublicUrl() + "/" + ANNOTATION_URL_PATH);
			return annotationUrl;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
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

	@Override
	public void RegionChanged(BpCoordRegion bpRegion) {
		long location = bpRegion.getMid();
		megaLocation.setText("" + (location/1000000));
		kiloLocation.setText("" + (location % 1000000) / 1000);
		unitLocation.setText("" + (location % 1000));
	}

	@Override
	public void vetoableChange(PropertyChangeEvent evt)
			throws PropertyVetoException {
		Object source = evt.getSource();
		if (source == megaLocation || source == kiloLocation || source == unitLocation) {
			plot.moveDataBpRegion(Long.parseLong(megaLocation.getText())*1000000 +
					Long.parseLong(megaLocation.getText())*1000 +
					Long.parseLong(megaLocation.getText()));
		}
	}
}
