package fi.csc.microarray.client.visualisation;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.methods.threed.DataPoint;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class AnnotateListPanel extends JPanel {

	private JList selectedList;
	private DefaultListModel selectedListModel;
	private JButton filterButton;
	private List<DataBean> datas = new ArrayList<DataBean>();

	private ClientApplication application = Session.getSession().getApplication();

	private JLabel countLabel;
	private String nameOfItems = "Genes";

	public AnnotateListPanel() {

		this.setLayout(new BorderLayout());
		selectedListModel = new DefaultListModel();
		selectedList = new JList(selectedListModel);

		// disabling selections changes text to gray, so we just make it look
		// like not selectable
		selectedList.setSelectionBackground(selectedList.getBackground());
		selectedList.setSelectionForeground(selectedList.getForeground());

		countLabel = new JLabel();

		filterButton = new JButton("Create dataset");
		filterButton.setToolTipText("Create new dataset from selected genes");
		filterButton.setEnabled(false);
		filterButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				VisualisationUtilities.filterBySelection(datas);
			}
		});

		this.add(countLabel, BorderLayout.NORTH);
		this.add(new JScrollPane(selectedList), BorderLayout.CENTER);

		JPanel buttonPanel = new JPanel(new BorderLayout());
		buttonPanel.add(filterButton, BorderLayout.SOUTH);
		this.add(buttonPanel, BorderLayout.SOUTH);
	}

	public AnnotateListPanel(String nameOfItems, boolean filterButtonVisible){
		this();
		this.nameOfItems  = nameOfItems;
		this.filterButton.setVisible(filterButtonVisible);
	}

	public void setSelectedListContentMultipleDatas(List<String> content, Map<DataBean, Set<Integer>> indexes, Object source, boolean showAnnotations, boolean dispatchEvent) {

		setData(indexes.keySet());

		selectedListModel.removeAllElements();
		setCount(content.size());
		filterButton.setEnabled(content.size() > 0);

		if (showAnnotations) {
			List<TableAnnotationProvider> annotationProviders = new LinkedList<TableAnnotationProvider>();
			try {
				for(DataBean data : indexes.keySet()){
					annotationProviders.add(new TableAnnotationProvider(data));
				}			
			} catch (MicroarrayException me) {
				throw new RuntimeException(me);
			}


			for (String row : content) {

				String viewName = "";
				for(TableAnnotationProvider annotation : annotationProviders){
					String annotationStr = annotation.getAnnotatedRowname(row);
					if(annotationStr != null && !viewName.contains(annotationStr)) {
						if(viewName.length() == 0){
							viewName += annotationStr;
						}else{
							viewName += ", " + annotationStr;
						}
					}	
				}

				selectedListModel.addElement(viewName);
			} 
		} else {
			for (String row : content) {
				selectedListModel.addElement(row);
			}
		}


		if (dispatchEvent) {
			for (DataBean data : indexes.keySet()) {
				application.getSelectionManager().getSelectionManager(data).setSelected(indexes.get(data), source);
			}
		}
	}

	/**
	 * TODO This method is only for 3d-scatter, as selection system is moving towards
	 * handling only row numbers. 3d-scatter should use this convention also and this 
	 * method could be removed after that. 
	 * 
	 * @param content
	 * @param source
	 * @param dispatchEvent
	 * @param data
	 *            is needed only if event is dispatched
	 */
	public void setSelectedListContentAsDataPoints(Collection<DataPoint> content, Object source, boolean dispatchEvent, DataBean data) {

		setData(data);

		TableAnnotationProvider annotationProvider;
		try {
			annotationProvider = new TableAnnotationProvider(data);

		} catch (MicroarrayException me) {
			throw new RuntimeException(me);
		}		

		selectedListModel.removeAllElements();
		setCount(content.size());
		filterButton.setEnabled(content.size() > 0);
		/*
		 * for(DataPoint row: content){
		 * selectedListModel.addElement(row.toString()); }
		 */

		int[] indexes = new int[content.size()];
		int i = 0;
		/*//Without gene symbols
		for (DataPoint row : content) {
			selectedListModel.addElement(row.toString());
			indexes[i++] = row.getIndex();
		}*/

		for (DataPoint row : content) {
			selectedListModel.addElement(annotationProvider.getAnnotatedRowname(row.toString()));
			indexes[i++] = row.getIndex();
		}


		if (dispatchEvent) {
			application.getSelectionManager().getSelectionManager(data).setSelection(indexes, source);
		}
	}

	/**
	 * @param rows
	 * @param source
	 * @param dispatchEvent
	 * @param data
	 *            is needed only if event is dispatched
	 */
	public void setSelectedRows(Set<Integer> rows, Object source, boolean dispatchEvent, DataBean data) {

		setData(data);
		TableAnnotationProvider annotationProvider;
		try {
			annotationProvider = new TableAnnotationProvider(data);

		} catch (MicroarrayException me) {
			throw new RuntimeException(me);
		}

		selectedListModel.removeAllElements();
		setCount(rows.size());
		filterButton.setEnabled(rows.size() > 0);


		//TODO getAnnotatedRowname should allow row index arguments, as it is used generally
		//to locate rows in chipster. After that finding these identifiers isn't necessary anymore
		Iterator<String> ids = null;
		try {
			ids = data.queryFeatures("/identifier").asStrings().iterator();
		} catch (MicroarrayException e) {
			//Finding identifiers shouldn't be necessary at all, see TODO couple rows upwards
			application.reportException(e);
		}

		for(int i = 0; ids.hasNext(); i++){
			String id = ids.next();

			if(rows.contains(i)){
				selectedListModel.addElement(annotationProvider.getAnnotatedRowname(id));
			}
		}

		if (dispatchEvent) {
			application.getSelectionManager().getSelectionManager(data).setSelected(
					rows, source);
		}
	}

	private void setCount(int count){
		countLabel.setText(count + " " + nameOfItems + " selected");
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
