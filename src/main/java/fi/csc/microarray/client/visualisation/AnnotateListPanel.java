package fi.csc.microarray.client.visualisation;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.methods.Scatterplot;
import fi.csc.microarray.client.visualisation.methods.threed.DataPoint;
import fi.csc.microarray.databeans.DataBean;

public class AnnotateListPanel extends JPanel {

	private JList selectedList;
	private DefaultListModel selectedListModel;
	private JButton annotateButton;
	private JButton filterButton;
	private List<DataBean> datas = new ArrayList<DataBean>();

	private ClientApplication application = Session.getSession().getApplication();

	private JLabel countLabel;

	public AnnotateListPanel() {

		this.setLayout(new BorderLayout());
		selectedListModel = new DefaultListModel();
		selectedList = new JList(selectedListModel);

		// disabling selections changes text to gray, so we just make it look
		// like it
		selectedList.setSelectionBackground(selectedList.getBackground());
		selectedList.setSelectionForeground(selectedList.getForeground());

		countLabel = new JLabel();

		annotateButton = new JButton("Annotate");
		annotateButton.setToolTipText("Create dataset and annotate with Bioconductor");
		annotateButton.setEnabled(false);
		annotateButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				VisualisationUtilities.annotateBySelection(datas);
			}
		});

		filterButton = new JButton("Create dataset");
		filterButton.setToolTipText("Create new dataset from selected rows");
		filterButton.setEnabled(false);
		filterButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				VisualisationUtilities.filterBySelection(datas);
			}
		});

		this.add(countLabel, BorderLayout.NORTH);
		this.add(new JScrollPane(selectedList), BorderLayout.CENTER);

		JPanel buttonPanel = new JPanel(new BorderLayout());
		buttonPanel.add(annotateButton, BorderLayout.NORTH);
		buttonPanel.add(filterButton, BorderLayout.SOUTH);
		this.add(buttonPanel, BorderLayout.SOUTH);
	}

	public void setSelectedListContentMultipleDatas(List<String> content, Map<DataBean, Set<Integer>> indexes, Object source, boolean dispatchEvent) {

		setData(indexes.keySet());

		selectedListModel.removeAllElements();
		countLabel.setText(content.size() + " Points selected");
		annotateButton.setEnabled(content.size() > 0);
		filterButton.setEnabled(content.size() > 0);

		for (String row : content) {
			selectedListModel.addElement(row.toString());
		}

		if (dispatchEvent) {
			for (DataBean data : indexes.keySet()) {
				application.getSelectionManager().getRowSelectionManager(data).setSelected(indexes.get(data), source);
			}
		}
	}

	/**
	 * @param content
	 * @param source
	 * @param dispatchEvent
	 * @param data
	 *            is needed only if event is dispatched
	 */
	public void setSelectedListContent(Collection<DataPoint> content, Object source, boolean dispatchEvent, DataBean data) {

		setData(data);
		selectedListModel.removeAllElements();
		countLabel.setText(content.size() + " Points selected");
		annotateButton.setEnabled(content.size() > 0);
		filterButton.setEnabled(content.size() > 0);
		/*
		 * for(DataPoint row: content){
		 * selectedListModel.addElement(row.toString()); }
		 */

		int[] indexes = new int[content.size()];
		int i = 0;
		for (DataPoint row : content) {
			selectedListModel.addElement(row.toString());
			indexes[i++] = row.getIndex();
		}

		if (dispatchEvent) {
			application.getSelectionManager().getRowSelectionManager(data).setSelected(indexes, source);
		}
	}

	/**
	 * @param content
	 * @param source
	 * @param dispatchEvent
	 * @param data
	 *            is needed only if event is dispatched
	 * @throws MicroarrayException
	 */
	public void setSelectedListContentAsDataItems(Collection<Scatterplot.DataItem2D> content, Object source, boolean dispatchEvent, DataBean data) {

		setData(data);
		TableAnnotationProvider annotationProvider;
		try {
			annotationProvider = new TableAnnotationProvider(data);
			
		} catch (MicroarrayException me) {
			throw new RuntimeException(me);
		}

		selectedListModel.removeAllElements();
		countLabel.setText(content.size() + " Points selected");
		annotateButton.setEnabled(content.size() > 0);
		filterButton.setEnabled(content.size() > 0);

		int[] indexes = new int[content.size()];
		int i = 0;
		for (Scatterplot.DataItem2D row : content) {
			selectedListModel.addElement(annotationProvider.getAnnotatedRowname(row.getName()));
			indexes[i++] = row.getRowIndex();
		}
		if (dispatchEvent) {
			application.getSelectionManager().getRowSelectionManager(data).setSelected(indexes, source);
		}
	}

	/**
	 * The datas are stored in a list to cope with visualisations for multiple
	 * datasets. This is a helper method to handle the list when there is only
	 * one dataset.
	 * 
	 * @param data
	 */
	private void setData(DataBean data) {
		datas.clear();
		datas.add(data);
	}

	private void setData(Set<DataBean> set) {
		datas.clear();
		for (DataBean data : set) {
			datas.add(data);
		}
	}
}
