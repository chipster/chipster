package fi.csc.microarray.databeans;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

/**
 * DataFolder is used to manage DataBean objects.
 * 
 * @see DataBean
 * @author Aleksi Kallio, hupponen
 * 
 */
public class DataFolder extends DataItemBase {

	public DataFolder(DataManager manager, String name) {
		this.manager = manager;
		this.name = name;
	}

	private List<DataItem> children = new LinkedList<DataItem>();
	private DataManager manager;

	public void addChild(DataItem child) {

		// was it already connected?
		boolean wasConnected = child.getParent() != null;

		// connect to this
		child.setParent(this);

		// add
		children.add(child);

		if (child instanceof DataBean) {
			doBackwardsCompatibleTypeTagInitialisation((DataBean) child);
		}

		// dispatch events if needed
		if (!wasConnected) {
			manager.dispatchEvent(new DataItemCreatedEvent(child));
		}
	}

	// TODO remove and rely on new type system
	private void doBackwardsCompatibleTypeTagInitialisation(DataBean data) {

		try {

			if (data.isContentTypeCompatitible("text/tab", "application/cel", "text/csv")) {
				data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
			}
			
			if (data.isContentTypeCompatitible("text/tab")) {
				BufferedReader in = null;
				try {
					in = new BufferedReader(new InputStreamReader(data.getContentByteStream()));
					String line = in.readLine();
					
					if (line != null) {						

						String[] split = line.split("\t");					

						// Some special tags for the tsv output of tool Annotate variants 
						if (split.length > 4 && 
								"SEQNAMES".equals(split[1]) && 
								"START".equals(split[2]) && 
								"END".equals(split[3])) {

							data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN);
							data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN);
							data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN);
						}

						// Cna data, first data column of rownames is missing in header
						if (split.length > 3 && 
								"chromosome".equals(split[0]) && 
								"start".equals(split[1]) && 
								"end".equals(split[2])) {

							data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN);
							data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN);
							data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN);
						}
					}
					
				} catch (IOException e) {
					throw new RuntimeException(e);
					
				} finally {
					IOUtils.closeIfPossible(in);
				}				
			}
			
			if (data.isContentTypeCompatitible("text/gtf")) {
				data.addTypeTag(MicroarrayModule.TypeTags.GTF_FILE);
				data.addTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
				data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
				data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_FOURTH_TABLE_COLUMN);				
				data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FIFTH_TABLE_COLUMN);				
				addTypeTagIfHashHeader(data);
			}
			
			if (data.isContentTypeCompatitible("text/vcf")) {
				data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
				addTypeTagIfDoubleHashHeader(data);
				data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
				data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN);
			}

			if (data.isContentTypeCompatitible("text/bed")) {
				data.addTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
				
				data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
				data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN);				
				data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN);

				// Check if it has title row
				BufferedReader in = null;
				try {
					in = new BufferedReader(new InputStreamReader(data.getContentByteStream()));
					String line = in.readLine();
					if (line != null && line.startsWith("track")) {
						data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_TITLE_ROW);
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				} finally {
					IOUtils.closeIfPossible(in);
				}
				
			}

			// mothur
			if (data.isContentTypeCompatitible("text/mothur-oligos")) {
				data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_OLIGOS);
			}
			if (data.isContentTypeCompatitible("text/mothur-names")) {
				data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_NAMES);
			}
			if (data.isContentTypeCompatitible("text/mothur-groups")) {
				data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_GROUPS);
			}
			
			
			// the rest is microarray specific
			if (!(Session.getSession().getPrimaryModule() instanceof MicroarrayModule)) {
				return;
			}

			Table chips = data.queryFeatures("/column/chip.*").asTable();

			// Tag the "main type"
			
			if (data.isContentTypeCompatitible("application/cel")) {
				data.addTypeTag(MicroarrayModule.TypeTags.RAW_AFFYMETRIX_EXPRESSION_VALUES);

			} else if (data.queryFeatures("/column/sample").exists() && !data.queryFeatures("/phenodata").exists()) {
				data.addTypeTag(MicroarrayModule.TypeTags.RAW_EXPRESSION_VALUES);

			} else if (chips != null && chips.getColumnCount() > 0) {
				data.addTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);
				
			} else if (data.isContentTypeCompatitible("application/bam")) {
				data.addTypeTag(MicroarrayModule.TypeTags.BAM_FILE);
				
			}  else if (data.isContentTypeCompatitible("chemical/x-fasta")) {
				data.addTypeTag(MicroarrayModule.TypeTags.FASTA_FILE);
			} 
			
			if (data.queryFeatures("/identifier").exists()) {
				data.addTypeTag(MicroarrayModule.TypeTags.GENENAMES);
			} 


			// Tag additional typing information
			if (data.queryFeatures("/phenodata").exists()) {
				data.addTypeTag(BasicModule.TypeTags.PHENODATA);
			}
				
			if (data.queryFeatures("/column/p.*").exists() && data.queryFeatures("/column/FC*").exists()) {
				data.addTypeTag(MicroarrayModule.TypeTags.SIGNIFICANT_EXPRESSION_FOLD_CHANGES);
			}

			boolean isChipwise = false;
			ParameterRecord pcaOn = data.getOperationRecord().getParameter("do.pca.on");
			if (pcaOn != null) {
				String pcaOnValue = pcaOn.getValue();
				if (pcaOnValue != null && pcaOnValue.equals("chips")) {
					isChipwise = true;
				}
			}
			if (data.getOperationRecord().getNameID().getID().equals("ordination-pca.R") && isChipwise) {
				data.addTypeTag(MicroarrayModule.TypeTags.EXPRESSION_PRIMARY_COMPONENTS_CHIPWISE);
			}

			if (chips != null && chips.getColumnNames().length > 1 && data.queryFeatures("/column/cluster").exists()) {
				data.addTypeTag(MicroarrayModule.TypeTags.CLUSTERED_EXPRESSION_VALUES);
			}
			
		    if (data.isContentTypeCompatitible("text/plain", "text/bed", "text/tab") 
		    		|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().contains(".bam-summary")) 
		    		|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().endsWith(".bam") || data.getName().endsWith(".sam"))
		    		|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().endsWith(".bai"))) {

		    	data.addTypeTag(MicroarrayModule.TypeTags.ORDERED_GENOMIC_ENTITIES);
		    }

		    if (data.queryFeatures("/clusters/som").exists()) {
		    	data.addTypeTag(MicroarrayModule.TypeTags.SOM_CLUSTERED_EXPRESSION_VALUES);
		    }
			
		} catch (MicroarrayException e) {
			throw new RuntimeException(e);
		}

	}

	public void removeChild(DataItem child) {
		// remove connections
		child.setParent(null);

		// remove
		children.remove(child);

		// dispatch events
		manager.dispatchEvent(new DataItemRemovedEvent(child));
	}

	public Iterable<DataItem> getChildren() {
		return children;
	}

	public int getChildCount() {
		return children.size();
	}

	public DataFolder getChildFolder(String name) {
		for (DataItem child : getChildren()) {
			if (child instanceof DataFolder && child.getName().equals(name)) {
				return (DataFolder) child;
			}
		}
		return null;
	}

	public String toString() {
		return getName();
	}

	private void addTypeTagIfHashHeader(DataBean data) {
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(data.getContentByteStream()));
			String line = in.readLine();
			if (line != null && line.startsWith("#")) {
				data.addTypeTag(MicroarrayModule.TypeTags.TABLE_WITH_HASH_HEADER);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeIfPossible(in);
		}
	}

	private void addTypeTagIfDoubleHashHeader(DataBean data) {
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(data.getContentByteStream()));
			String line = in.readLine();
			if (line != null && line.startsWith("##")) {
				data.addTypeTag(MicroarrayModule.TypeTags.TABLE_WITH_DOUBLE_HASH_HEADER);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeIfPossible(in);
		}
	}


}
