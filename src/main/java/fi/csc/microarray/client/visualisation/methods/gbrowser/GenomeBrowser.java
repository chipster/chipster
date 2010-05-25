package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.PlainDocument;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Row;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.util.IOUtils;

/**
 * @author Petri Klemelä, Aleksi Kallio
 */
public class GenomeBrowser extends Visualisation implements ActionListener, RegionListener, FocusListener {

	private static final String[] CHROMOSOMES = new String[] {
		"1",
		"2",
		"3",
		"4",
		"5",
		"6",
		"7",
		"8",
		"9",
		"10",
		"11",
		"12",
		"13",
		"14",
		"15",
		"16",
		"17",
		"18",
		"19",
		"20",
		"21",
		"22",
		"X",
		"Y",
	};
	
	private static final long[] CHROMOSOME_SIZES = new long[] {
		247199719L,	
		242751149L,
		199446827L,
		191263063L,
		180837866L,
		170896993L,
		158821424L,
		146274826L, 	
		140442298L,
		135374737L,
		134452384L,
		132289534L, 	
		114127980L,
		106360585L,
		100338915L,
		88822254L, 	
		78654742L,
		76117153L,
		63806651L,
		62435965L,
		46944323L,
		49528953L,
		154913754L,
		57741652L, 	
	};
	private static final String ANNOTATION_URL_PATH = "annotations";
	private static final String CONTENTS_FILE = "contents.txt";

	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";

	private static enum TrackType {
		CYTOBANDS, GENES, TRANSCRIPTS, TREATMENT_READS, CONTROL_READS, PEAKS, REFERENCE
	}

	private static class Track {
		
		TrackType type;
		JCheckBox checkBox; 
		String name;
		DataBean userData;

		public Track(String name, TrackType type) {
			this.name = name;
			this.type = type;
			
		}

		public Track(String name, TrackType type, DataBean userData) {
			this(name, type);
			this.userData = userData;
		}
	}

	private final ClientApplication application = Session.getSession().getApplication();

	private List<DataBean> datas;
	private List<Track> tracks = new LinkedList<Track>();

	private GenomePlot plot;

	private JPanel paramPanel;
	private JPanel settingsPanel = new JPanel();
	private JPanel plotPanel = new JPanel(new CardLayout());

	private JButton gotoButton = new JButton("Go to location");
	private JButton drawButton = new JButton("Draw");

	private JTextField megaLocation = new JTextField(4);
	private JTextField kiloLocation = new JTextField(4);
	private JTextField unitLocation = new JTextField(4);
	private JComboBox chrBox = new JComboBox();
	private JComboBox genomeBox = new JComboBox();
	// private JRadioButton horizView;
	// private JRadioButton circularView;
	private GridBagConstraints settingsGridBagConstraints;
	private List<Row> contents;

	private String localAnnotationPath;

	private URL annotationUrl;


	public GenomeBrowser(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public JPanel getParameterPanel() {

		// FIXME should the following check be enabled?
		if (paramPanel == null /* || data != application.getSelectionManager().getSelectedDataBean() */) {

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

	private void createAvailableTracks() {
		
		String genome = (String)genomeBox.getSelectedItem();
		
		// list available track types for the genome
		for (Row row : contents) {
			if (genome.equals(row.version)) {
				TrackType type;
				if (row.content.contains("Genes")) {
					type = TrackType.GENES;					
				} else if (row.content.contains("Transcripts")) {
					type = TrackType.TRANSCRIPTS;
				} else if (row.content.contains("Cytobands")) {
					type = TrackType.CYTOBANDS;
				} else if (row.content.contains("Reference")) {
					continue; // track not directly supported, skip
				} else {
					continue; // track not supported, skip
				}

				tracks.add(new Track(row.content, type));
			}
		}
		
		List<TrackType> interpretations = interpretUserDatas(this.datas);

		for (int i = 0; i < interpretations.size(); i++) {
			TrackType interpretation = interpretations.get(i);
			tracks.add(new Track(datas.get(i).getName(), interpretation, datas.get(i)));
		}

		// list available track types for the genome
		for (Track track : tracks) {
			this.settingsGridBagConstraints.gridy++;
			JCheckBox box = new JCheckBox(track.name, true);
			settingsPanel.add(box, this.settingsGridBagConstraints);
			track.checkBox = box;
		}
		
		GridBagConstraints c = this.settingsGridBagConstraints;
		c.gridy++;
		settingsPanel.add(drawButton, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		settingsPanel.add(new JPanel(), c);

	}

	public JPanel createSettingsPanel() {

		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

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

		for (JTextField field : new JTextField[] {megaLocation, kiloLocation, unitLocation}) {
			PlainDocument fieldContents = new PlainDocument() {
				@Override
				public void insertString(int offs, String str, AttributeSet a) throws BadLocationException {
					if (str != null && str.length() > 3) {
						return; // was too long
					}
					super.insertString(offs, str, a);
				}
			};			
			field.setDocument(fieldContents);
			field.addFocusListener(this);
		}

		settingsPanel.add(new JLabel("Location"), c);

		c.anchor = GridBagConstraints.SOUTH;
		JLabel megaLabel = new JLabel("M");
		JLabel kiloLabel = new JLabel("k");
		c.gridy++;
		c.gridwidth = 1;
		c.insets.set(5, 10, 5, 0);
		c.weightx = 1.0;
		settingsPanel.add(megaLocation, c);
		c.gridx++;
		c.insets.set(5, 0, 5, 0);
		c.weightx = 0.0;
		settingsPanel.add(megaLabel, c);
		c.gridx++;
		c.weightx = 1.0;
		settingsPanel.add(kiloLocation, c);
		c.gridx++;
		c.weightx = 0.0;
		settingsPanel.add(kiloLabel, c);
		c.gridx++;
		c.insets.set(5, 0, 5, 10);
		c.weightx = 1.0;
		settingsPanel.add(unitLocation, c);

		c.gridx = 0;
		c.gridwidth = 5;		
		c.gridy++;
		c.insets.set(5, 10, 5, 10);
		settingsPanel.add(gotoButton , c);
		gotoButton.addActionListener(this);
		gotoButton.setEnabled(false);

		InputStream contentsStream = null;
		
		try {
			// parse what annotations we have available
			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
    		String configuredLocalPath = configuration.getString("client", "local-annotation-path");
			if (configuredLocalPath.trim().isEmpty()) {
				this.localAnnotationPath = null;
				this.annotationUrl = fetchAnnotationUrl();
				contentsStream = new URL(annotationUrl + "/" + CONTENTS_FILE).openStream();
			} else {
				this.localAnnotationPath = configuredLocalPath;
				this.annotationUrl = null;
				contentsStream = new FileInputStream(localAnnotationPath + File.separator + CONTENTS_FILE);
			}
			
			AnnotationContents annotationContentFile = new AnnotationContents();
			annotationContentFile.parseFrom(contentsStream);
			this.contents = annotationContentFile.getRows();
		
			// read genome name and version for each annotation file
			LinkedHashSet<String> genomes = annotationContentFile.getGenomes();
			c.gridy++;
			settingsPanel.add(new JLabel("Genome"), c);
			c.gridy++;
			for (String genome : genomes) {
				genomeBox.addItem(genome);
			}
			settingsPanel.add(genomeBox, c);
		
			// list available chromosomes
			// FIXME These should be read from user data file
			for (String chromosome : CHROMOSOMES) {
				chrBox.addItem(chromosome);
			}
		
			c.gridy++;
			settingsPanel.add(new JLabel("Chromosome"), c);
			c.gridy++;
			settingsPanel.add(chrBox, c);
			
		} catch (IOException e) {
			application.reportException(e);
		
		} finally {
			IOUtils.closeIfPossible(contentsStream);
		}

		// horizView = new JRadioButton("Horizontal");
		// horizView.setSelected(true);
		// circularView = new JRadioButton("Circular");
		//
		// views = new ButtonGroup();
		// views.add(horizView);
		// views.add(circularView);
		//
		// c.gridy++;
		// settingsPanel.add(new JLabel("View mode"), c);
		// c.gridy++;
		// settingsPanel.add(horizView, c);
		// c.gridy++;
		// settingsPanel.add(circularView, c);

		this.settingsGridBagConstraints = c;
		

		
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
		} else if (source == gotoButton) {
			gotoButton.setEnabled(false);
			locationChanged();
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		return getVisualisation(Arrays.asList(new DataBean[] { data }));
	}

	@Override
	public JComponent getVisualisation(java.util.List<DataBean> datas) throws Exception {
		this.datas = datas;
		createAvailableTracks(); // we can create tracks now that we know the data
		
		// create panel with card layout and put message panel there
		JPanel waitPanel = new JPanel();
		waitPanel.add(new JLabel("Please select parameters"));
		plotPanel.add(waitPanel, WAITPANEL);

		return plotPanel;
	}

	private void showVisualisation() {

		try {
			// create the plot
			String genome = (String) genomeBox.getSelectedItem();
			this.plot = new GenomePlot(true);

			// add selected tracks
			LinkedList<DataSource> treatments = new LinkedList<DataSource>();
			LinkedList<DataSource> controls = new LinkedList<DataSource>();
			LinkedList<DataSource> peaks = new LinkedList<DataSource>();
			
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					File file = track.userData == null ? null : Session.getSession().getDataManager().getLocalFile(track.userData);
					switch (track.type) {
					case CYTOBANDS:
						TrackFactory.addCytobandTracks(plot, createAnnotationDataSource("Homo_sapiens.GRCh37.57_karyotype.tsv")); // using always the
						break;
					case GENES:
						TrackFactory.addGeneTracks(plot, createAnnotationDataSource("Homo_sapiens." + genome + "_genes.tsv"));
						break;
					case REFERENCE:
						// integrated into peaks
						break;
					case TRANSCRIPTS:
						TrackFactory.addTranscriptTracks(plot, createAnnotationDataSource("Homo_sapiens." + genome + "_transcripts.tsv"));
						break;
					case PEAKS:
						peaks.add(new DataSource(file));
						break;
					case TREATMENT_READS:
						treatments.add(new DataSource(file));
						break;
					case CONTROL_READS:
						controls.add(new DataSource(file));
						break;
					}
				}
			}

			// add user tracks and ruler
			if (!peaks.isEmpty()) {
				TrackFactory.addPeakTracks(plot, peaks);
			}
			if (!treatments.isEmpty() || !controls.isEmpty()) {
				TrackFactory.addReadTracks(plot, treatments, controls, createAnnotationDataSource("Homo_sapiens." + genome + "_seq.tsv"));
			}
			TrackFactory.addRulerTrack(plot);

			// initialise the plot
			plot.addDataRegionListener(this);
			plot.start((String)chrBox.getSelectedItem(), (double)CHROMOSOME_SIZES[chrBox.getSelectedIndex()]);
			
			// wrap it in a panel
			ChartPanel chartPanel =  new NonScalableChartPanel(new JFreeChart(plot));
			plot.chartPanel = chartPanel;
			chartPanel.setCursor(new Cursor(Cursor.HAND_CURSOR));
			
			// add mouse listeners
			for (View view : plot.getViews()) {
				chartPanel.addMouseListener(view);
				chartPanel.addMouseMotionListener(view);
				chartPanel.addMouseWheelListener(view);
			}

			// put panel on top of card layout
			if (plotPanel.getComponentCount() == 2) {
				plotPanel.remove(1);
			}
			plotPanel.add(chartPanel, PLOTPANEL);
			CardLayout cl = (CardLayout) (plotPanel.getLayout());
			cl.show(plotPanel, PLOTPANEL);
			
		} catch (Exception e) {
			application.reportException(e);
		}
	}

	private DataSource createAnnotationDataSource(String file) throws FileNotFoundException, MalformedURLException {
		if (this.annotationUrl != null) {
			return new DataSource(this.annotationUrl, file);
		} else {
			return new DataSource(new File(this.localAnnotationPath), file);
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
	public boolean canVisualise(DataBean data) throws MicroarrayException {
		return canVisualise(Arrays.asList(new DataBean[] { data }));
	}

	@Override
	public boolean canVisualise(java.util.List<DataBean> datas) throws MicroarrayException {
		return interpretUserDatas(datas) != null;
	}

	public class ObjVariable extends Variable {

		public Object obj;

		public ObjVariable(Object obj) {
			super(null, null);
			this.obj = obj;
		}
	}

	@Override
	public void regionChanged(BpCoordRegion bpRegion) {
		long location = bpRegion.getMid();
		megaLocation.setText("" + (location / 1000000));
		kiloLocation.setText("" + (location % 1000000) / 1000);
		unitLocation.setText("" + (location % 1000));
		gotoButton.setEnabled(false);
	}

	private List<TrackType> interpretUserDatas(List<DataBean> datas) {
		LinkedList<TrackType> interpretations = new LinkedList<TrackType>();

		// try to find interpretation for all selected datas
		for (DataBean data : datas) {

			if (data.isContentTypeCompatitible("text/plain")) {
				// reads
				if (data.getName().contains("control")) {
					interpretations.add(TrackType.CONTROL_READS);
				} else {
					interpretations.add(TrackType.TREATMENT_READS);
				}

			} else if (data.isContentTypeCompatitible("text/bed")) {
				// peaks
				interpretations.add(TrackType.PEAKS);

			} else {
				// cannot interpret, visualisation not available for this selection
				return null;
			}
		}

		return interpretations;
	}

	@Override
	public boolean isForSingleData() {
		return true;
	}

	@Override
	public boolean isForMultipleDatas() {
		return true;
	}

	private void locationChanged() {
		plot.moveDataBpRegion(Long.parseLong(megaLocation.getText()) * 1000000 + Long.parseLong(kiloLocation.getText()) * 1000 + Long.parseLong(unitLocation.getText()));
	}

	@Override
	public void focusGained(FocusEvent e) {
		gotoButton.setEnabled(true);		
	}

	@Override
	public void focusLost(FocusEvent e) {
		// skip		
	}
}
