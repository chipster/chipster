package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TranscriptTrack extends Track{

	private Map<String, Gene> genes = 
		new TreeMap<String, Gene>();

	List<Integer> occupiedSpace = new ArrayList<Integer>();

	private Color color;

	private int RESOLUTION = 512;
		
	public enum PartColor { CDS (new Color(64, 192, 64)), UTR(new Color(192, 64, 64)), START_CODON(Color.gray);
		public Color c;

		PartColor(Color c){
			this.c = c;
		}
	}

	public TranscriptTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			FileParser inputParser, Color color, long maxBpLength)  {

		super(view, file, handler, inputParser);		
		this.color = color;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		if(genes != null){

			List<Gene> sortedGenes = new ArrayList<Gene>(genes.values());
			Collections.sort(sortedGenes);

			for(Gene gene : sortedGenes){

				if(!gene.region.intercepts(getView().getBpRegion())){
					
					genes.remove(gene.id);
					continue;
				}

				Rectangle rect = new Rectangle();

				rect.x = getView().bpToTrack(gene.region.start);
				rect.width = getView().bpToTrack(gene.region.end) - rect.x;

				int i = 0;

				while(occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x ){
					i++;
				}

				int end = rect.x + rect.width;

				if(occupiedSpace.size() > i){
					occupiedSpace.set(i, end);
				} else {
					occupiedSpace.add(end);
				}

				rect.y = (int)(getView().getTrackHeight() - ((i + 1) * (14)));
				rect.height = 2;

				rect.y += 1;
				drawables.add(new RectDrawable(rect, Color.darkGray, null));
				rect.y -= 1;
				
				rect.height = 4;
				
				//Draw arrow
				if(gene.first().values.get(ColumnType.STRAND) == Strand.REVERSED){
					drawables.addAll(
							getArrowDrawables(rect.x, rect.y, -rect.height, rect.height));
				} else {
					drawables.addAll(
							getArrowDrawables(rect.x + rect.width, rect.y, rect.height, rect.height));
				}
				
				String geneId = ((String)gene.first().values.get(ColumnType.DESCRIPTION));
				
				if(rect.width > geneId.length() * 5){
					drawables.add(new TextDrawable(
							rect.x, rect.y - 1, geneId, 
							Color.DARK_GRAY));
				}

				List<Drawable> geneDrawables = new ArrayList<Drawable>();

				for(RegionContent part : gene){

					if(part.values == null){
						drawables.add(createDrawable(part.region.start, part.region.end, color));
					} else {


						String value = ((String) part.values.get(ColumnType.VALUE)).trim();
						Color c;

						if(value.equals("CDS")){						
							c = PartColor.CDS.c;
						} else if ( value.equals("exon")){
							c = PartColor.UTR.c;
						} else if ( value.equals("start_codon")){
							c = PartColor.START_CODON.c;
						} else if ( value.equals("stop_codon")){
							
							//TODO Check how this should be visualised
							c = PartColor.UTR.c;
						}
						else {
							System.out.println("Gene description not recognised: " + value);
							c = Color.blue;
						}												

						rect.x = getView().bpToTrack(part.region.start);						
						rect.width = getView().bpToTrack(part.region.end) - rect.x;
						rect.height = 4;

						geneDrawables.add(new RectDrawable(rect, c, null));												
						//drawables.add(new RectDrawable(rect, c, c));					
					}							
				}

				Collections.sort(geneDrawables, new Comparator<Drawable>(){
					public int compare(Drawable one, Drawable other) {

						if(one.color.equals(PartColor.CDS.c) && 
								other.color.equals(PartColor.UTR.c)){
							return 1;
						} else if ( one.color.equals(PartColor.UTR.c) && 
								other.color.equals(PartColor.CDS.c)){
							return -1;							
						} else {
							return 0;
						}																																	
					}
				});
				
				drawables.addAll(geneDrawables);
			}
		}
		return drawables;
	}

	private Drawable createDrawable(BpCoord startBp, BpCoord endBp, Color c){
		return createDrawable(startBp, endBp, 5, c);
	}

	private Drawable createDrawable(BpCoord startBp, BpCoord endBp, int height, Color c){
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;

		int i = 0;

		while(occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x + 1){
			i++;
		}

		int end = rect.x + rect.width;

		if(occupiedSpace.size() > i){
			occupiedSpace.set(i, end);
		} else {
			occupiedSpace.add(end);
		}

		rect.y = (int)(getView().getTrackHeight() - ((i + 1) * (height + 2)));
		rect.height = height;

		return new RectDrawable(rect, c, null);
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {		

		//Genes and transcripts are ordered in the file, but to here they come in any order
		//That's why we have to put them to Gene objects to sort them again
		//Sorting is needed to draw partly overlapping genes in the same order every time
		if(!areaResult.status.concise && areaResult.content.values.get(ColumnType.STRAND) == getStrand()){
			
			Map<ColumnType, Object> values = areaResult.content.values;
			String id = (String)values.get(ColumnType.PARENT_ID);
			
			if(!genes.containsKey(id)){
				genes.put(id, new Gene(new BpCoordRegion(
						(Long)values.get(ColumnType.PARENT_BP_START),
						(Long)values.get(ColumnType.PARENT_BP_END),
						(Chromosome)values.get(ColumnType.CHROMOSOME)), id));
			}

			genes.get(id).add(areaResult.content);

			getView().redraw();
		}
	}

	private boolean wasLastConsied = true;

	private long maxBpLength;

	@Override
	public void updateData(){

		if(wasLastConsied != isConcised()){
			genes.clear();
			wasLastConsied = isConcised();
		}
		super.updateData();
	}

	public int getMaxHeight(){
		if(getView().getBpRegion().getLength() <= maxBpLength){
			return super.getMaxHeight();
		} else {
			return 0;
		}
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] { 
				ColumnType.CHROMOSOME,
				ColumnType.PARENT_BP_START,		
				ColumnType.PARENT_BP_END,
				ColumnType.STRAND,
				ColumnType.DESCRIPTION,
				ColumnType.VALUE,							
				ColumnType.PARENT_ID,
				ColumnType.PARENT_PART
		});
	}

	@Override
	public boolean isConcised() {
		return false;
		//return getView().getBpRegion().getLength() > 1*1024*1024;
		//return reads.size() > RESOLUTION;
	}
}
