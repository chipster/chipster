package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.selection.RowChoiceEvent;
import fi.csc.microarray.client.selection.RowSelectionManager;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.VenndiPlot.AREAS;
import fi.csc.microarray.databeans.DataBean;

public class VennDiagram extends Visualisation implements PropertyChangeListener {

	private static final String IDENTIFIER_COLUMN = "/identifier";

	public VennDiagram(VisualisationFrame frame) {
		super(frame);
	}

	private VenndiPlot plot;
	private ChartPanel chartPanel;

	private JPanel paramPanel;
	private AnnotateListPanel list;

	@Override
	public JPanel getParameterPanel() {
		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			paramPanel.setLayout(new BorderLayout());

			list = new AnnotateListPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Selected", list);

			paramPanel.add(tabPane, BorderLayout.CENTER);
		}
		return paramPanel;
	}

	@Override
	public JComponent getVisualisation(List<DataBean> datas) throws Exception {

		if (datas.size() < 2 || datas.size() > 3) {
			throw new IllegalArgumentException("Venn Diagram can be used only with two or three datasets");
		}

		Map<String, Integer> A = new HashMap<String, Integer>();
		Map<String, Integer> B = new HashMap<String, Integer>();
		Map<String, Integer> C = new HashMap<String, Integer>();

		int i = 0;
		for (String name : datas.get(0).queryFeatures(IDENTIFIER_COLUMN).asStrings()) {
			A.put(name, i++);

		}

		i = 0;
		for (String name : datas.get(1).queryFeatures(IDENTIFIER_COLUMN).asStrings()) {
			B.put(name, i++);
		}

		i = 0;
		if (datas.size() == 3) {
			for (String name : datas.get(2).queryFeatures(IDENTIFIER_COLUMN).asStrings()) {
				C.put(name, i++);
			}
		}

		// Every AREA is mapped to the set of String containing the identifiers
		// of each area.
		Map<AREAS, Set<String>> sortedIds = new HashMap<AREAS, Set<String>>(VenndiPlot.AREAS.values().length);

		for (AREAS area : AREAS.values()) {
			sortedIds.put(area, new HashSet<String>());
		}

		sortedIds.get(AREAS.ABC).addAll(A.keySet());
		sortedIds.get(AREAS.ABC).retainAll(B.keySet());
		sortedIds.get(AREAS.ABC).retainAll(C.keySet());

		sortedIds.get(AREAS.AB).addAll(A.keySet());
		sortedIds.get(AREAS.AB).retainAll(B.keySet());
		sortedIds.get(AREAS.AB).removeAll(C.keySet());

		sortedIds.get(AREAS.BC).addAll(B.keySet());
		sortedIds.get(AREAS.BC).retainAll(C.keySet());
		sortedIds.get(AREAS.BC).removeAll(A.keySet());

		sortedIds.get(AREAS.AC).addAll(A.keySet());
		sortedIds.get(AREAS.AC).retainAll(C.keySet());
		sortedIds.get(AREAS.AC).removeAll(B.keySet());

		sortedIds.get(AREAS.A).addAll(A.keySet());
		sortedIds.get(AREAS.A).removeAll(B.keySet());
		sortedIds.get(AREAS.A).removeAll(C.keySet());

		sortedIds.get(AREAS.B).addAll(B.keySet());
		sortedIds.get(AREAS.B).removeAll(A.keySet());
		sortedIds.get(AREAS.B).removeAll(C.keySet());

		sortedIds.get(AREAS.C).addAll(C.keySet());
		sortedIds.get(AREAS.C).removeAll(A.keySet());
		sortedIds.get(AREAS.C).removeAll(B.keySet());

		String[][] idTable = new String[AREAS.values().length][];

		for (AREAS area : AREAS.values()) {
			idTable[area.ordinal()] = sortedIds.get(area).toArray(new String[0]);
		}

		Map<DataBean, Map<String, Integer>> indexMaps = new HashMap<DataBean, Map<String, Integer>>();

		indexMaps.put(datas.get(0), A);
		indexMaps.put(datas.get(1), B);
		if (datas.size() > 2) {
			indexMaps.put(datas.get(2), C);
		}

		plot = new VenndiPlot(idTable, datas, indexMaps, this);

		chartPanel = makePanel(new JFreeChart("Venn-diagram", plot));

		chartPanel.addChartMouseListener(plot);

		this.updateSelectionsFromApplication(false);

		application.addPropertyChangeListener(this);

		return chartPanel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return false;
	}

	@Override
	public boolean canVisualise(List<DataBean> beans) throws MicroarrayException {
		if (beans.size() < 2 || beans.size() > 3) {
			return false;
		}

		for (DataBean data : beans) {
			boolean isTabular = VisualisationMethod.SPREADSHEET.getHeadlessVisualiser().canVisualise(data);

			if (!(isTabular && data.queryFeatures(IDENTIFIER_COLUMN).exists())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isForSingleData() {
		return false;
	}

	@Override
	public boolean isForMultipleDatas() {
		return true;
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof RowChoiceEvent && evt.getSource() != this) {

			updateSelectionsFromApplication(false);
		}
	}

	private void updateSelectionsFromApplication(boolean dispatchEvent) {

		List<String> selected = new ArrayList<String>();
		try {

			for (DataBean data : plot.getDataset().getDataBeans()) {
				RowSelectionManager manager = application.getSelectionManager().getRowSelectionManager(data);
				selected.addAll(manager.getSelectedIdentifiers());
			}

		} catch (MicroarrayException e) {
			application.reportException(new MicroarrayException("Unable to get selected identifiers", e));
		}

		plot.setSelected(selected);
		plot.setSelectedForList(false);
		// To show new selections
		chartPanel.repaint();
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		throw new MicroarrayException("Venn Diagram can be used only with two or three datasets");
	}

	/**
	 * AnnotateListPanel.setSelectedListContentMultipleDatas can't be seen from
	 * the VenndiPlot, so we need to forward this information.
	 * 
	 * @param ids
	 * @param indexes
	 * @param venndiPlot
	 * @param dispatchEvent
	 */
	public void setSelectedListContent(List<String> ids, Map<DataBean, Set<Integer>> indexes, VenndiPlot venndiPlot, boolean dispatchEvent) {

		list.setSelectedListContentMultipleDatas(ids, indexes, this, dispatchEvent);
	}
}
