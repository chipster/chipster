package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

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

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class GeneTrack extends Track{

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

	public GeneTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
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

				//Zetas to make all chromosomes to be "smaller" in alphatetical order of Strings
				BpCoord minBp = new BpCoord(Long.MAX_VALUE, new Chromosome(Chromosome.MAX_VALUE));
				BpCoord maxBp = new BpCoord(Long.MIN_VALUE, new Chromosome(Chromosome.MIN_VALUE));
				String id = null;

				for(RegionContent part : gene){

					Object valueObj = part.values.get(ColumnType.DESCRIPTION);

					minBp = minBp.min(part.region.start);
					maxBp = maxBp.max(part.region.end);
					id  = (String)part.values.get(ColumnType.ID);
				}

				if(!(new BpCoordRegion(minBp, maxBp)).intercepts(getView().getBpRegion())){
					genes.remove(id);
					continue;
				}

				Rectangle rect = new Rectangle();

				rect.x = getView().bpToTrack(minBp);
				rect.width = getView().bpToTrack(maxBp) - rect.x;

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

				rect.y = (int)(getView().getTrackHeight() - ((i + 1) * (20 + 2)));
				rect.height = 2;

				rect.y += 4;
				drawables.add(new RectDrawable(rect, Color.darkGray, null));
				rect.y -= 4;
				
				rect.height = 10;
				
				//Draw arrow
				if(gene.first().values.get(ColumnType.STRAND) == Strand.REVERSED){
					drawables.addAll(
							getArrowDrawables(rect.x, rect.y, -rect.height, rect.height));
				} else {
					drawables.addAll(
							getArrowDrawables(rect.x + rect.width, rect.y, rect.height, rect.height));
				}
				
				String geneId = ((String)gene.first().values.get(ColumnType.ID)).split("\"")[1];
				
				if(rect.width > geneId.length() * 7){
					drawables.add(new TextDrawable(
							rect.x, rect.y, geneId, 
							Color.DARK_GRAY));
				}

				List<Drawable> geneDrawables = new ArrayList<Drawable>();

				for(RegionContent part : gene){

					if(part.values == null){
						drawables.add(createDrawable(part.region.start, part.region.end, color));
					} else {


						String value = ((String) part.values.get(ColumnType.DESCRIPTION)).trim();
						Color c;
						int height = 0;

						if(value.equals("CDS")){						
							c = PartColor.CDS.c;
						} else if ( value.equals("exon")){
							c = PartColor.UTR.c;
						} else if ( value.equals("start_codon")){
							c = PartColor.START_CODON.c;
						} else {
							System.out.println("Gene description not recognised: " + value);
							c = Color.blue;
						}												

						rect.x = getView().bpToTrack(part.region.start);						
						rect.width = getView().bpToTrack(part.region.end) - rect.x;
						rect.height = 10;

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

		String id = (String)areaResult.content.values.get(ColumnType.PARENT_ID);

		if(id != null && areaResult.content.values.get(ColumnType.STRAND) == getStrand()){
			if(!genes.containsKey(id)){
				genes.put(id, new Gene());
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
		return Arrays.asList(new ColumnType[] { ColumnType.DESCRIPTION, ColumnType.STRAND, ColumnType.ID}); 
	}

	@Override
	public boolean isConcised() {
		return false;
		//return getView().getBpRegion().getLength() > 1*1024*1024;
		//return reads.size() > RESOLUTION;
	}
}
