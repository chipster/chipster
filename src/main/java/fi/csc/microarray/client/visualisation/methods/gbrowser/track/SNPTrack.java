package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * single nucleotide polymorphism track
 * there are two ways of showing data:
 * 1. by nucleotides color
 * 2. by consequence to transcript
 * 
 * in the consequence to transcript way the color is get by simple prioritization
 * mechanizm. There is array of colors, which are ordered by importance and
 * consequences enumeration, which are also ordered by importance, and particular
 * color refer to particular consequence name.
 * 
 * If importance sequence is changed, so the colors array must be changed accordingly
 * 
 * @author zukauska
 *
 */

public class SNPTrack extends Track {
	
	// prioritize (1. stop, 2. frameshift, synonymous, non_synonymous)
	
	private static int width = 8;
	Collection<RegionContent> values = new LinkedList<RegionContent>();
	Long lastPosition;
	Color a = new Color(64, 192, 64, 128);
	Color c = new Color(64, 64, 192, 128);
	Color g = new Color(128, 128, 128, 128);
	Color t = new Color(192, 64, 64, 128);
	
	//ordered by prioritiazation
	private Color[] colors = new Color[] {
			new Color(255, 0, 0, 128),//red
			new Color(255, 0, 0, 128),//red
			new Color(255, 105, 180, 128),//pink
			new Color(255, 215, 0, 128),//gold
			new Color(50, 205, 50, 128),//lemon green
			new Color(190, 190, 190, 128),//gray
	};
	
	private int lastColorIndex = colors.length-1;
	
	boolean changeView = false;
	
	//ordered according colors
	enum ConsequenceName {
		STOP_GAINED, STOP_LOST, FRAMESHIFT_CODING, SYNONYMOUS_CODING,
		NON_SYNONYMOUS_CODING, UPSTREAM, INTRONIC, INTERGENETIC,  
		COMPLEX_INDEL, PARTIAL_CODON, REGULATORY_REGION, WITHIN_MATURE_mIRNA,
		PRIME5_UTR, PRIME3_UTR, NONE,  
	}

	public SNPTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler) {
		super(view, file, handler);
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		Iterator<RegionContent> iter = values.iterator();
		
		if (values != null) {
			while (iter.hasNext()) {
				RegionContent value = iter.next();
				
				// remove those that are not in this view
	            if (!value.region.intercepts(getView().getBpRegion())) {
	                iter.remove();
	                continue;
	            }
	            
	            if (lastPosition == null) {
					lastPosition = (long)getView().bpToTrack(value.region.start);
				}
	            String allele = (String)value.values.get(ColumnType.ALLELE);
	            long position = getView().bpToTrack(value.region.start);
	            if (changeView) {
	            	String consequence = (String)value.values.get(ColumnType.CONSEQUENCE_TO_TRANSCRIPT);
	            	
	            	//if we have multiple lines with same location
	            	if (lastPosition == position) {
	            		Iterator<Drawable> i = drawables.iterator();
	            		int j = 0;
	            		Drawable lastDrawable = null;
	            		
	            		//get the last drawable
	            		while (i.hasNext()) {
	            			lastDrawable = i.next();
	            			j++;
	            		}
	            		if (lastDrawable == null) {
	            			
	            		} else {
	            			int old = getColorNumber(lastDrawable.color);
	            			int now = ConsequenceName.valueOf(consequence).ordinal();
	            			
	            			//if importance is less than the old one, then just skip it
	            			if (now >= old) {
	            			} else {
	            				//change the color
	            				drawables.remove(lastDrawable);
	            				drawables.add(new RectDrawable((int)position, 1, width,
	            						getHeight(), colors[now], colors[now]));
	            			}
	            		}
	            	} else {
	            		drawables.add(new RectDrawable((int)position, 1, 
            				width, getHeight(), colors[getColorNumber(consequence)],
            				colors[getColorNumber(consequence)]));
	            	}
	            } else {
	            	//frameshift
            		if (allele.matches("[A-Z]/[A-Z]/-") || 
	            			allele.matches("[A-Z]/-") || allele.matches("[A-Z]/-/[A-Z]")) {
	            		drawables.add(new RectDrawable((int)position, 1, width, getHeight(), 
	            				colors[getColorNumber("FRAMESHIFT_CODING")], 
	            				colors[getColorNumber("FRAMESHIFT_CODING")]));
	            		
	            	} else if (allele.matches("[A-Z]/A")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), a, a));
		            } else if (allele.matches("[A-Z]/C")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), c, c));
		            } else if (allele.matches("[A-Z]/G")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), g, g));
		            } else if (allele.matches("[A-Z]/T")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), t, t));
		            }
	            }
	            lastPosition = position;
			}
		}
		
		return drawables;
	}
	
	private int getColorNumber(Color co) {
		int index = -1;
		for (Color c : colors) {
			index++;
			if (c.equals(co)) {
				break;
			}
		}
		if(index == -1) {
			index = lastColorIndex;
		}
		return index;
	}
	
	private int getColorNumber(String cn) {
		int index = ConsequenceName.valueOf(cn).ordinal();
		if (index > lastColorIndex) {
			index = lastColorIndex;
		}
		return index;
	}
	
	@Override
	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		if (areaResult.content.values.get(ColumnType.STRAND) == getStrand()) {
			values.add(areaResult.content);
		}
	}
	
	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new
        HashMap<DataSource, Set<ColumnType>>();
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.POSITION,
                ColumnType.STRAND,
                ColumnType.CONSEQUENCE_TO_TRANSCRIPT,
                ColumnType.ALLELE})));
		return datas;
	}
	
	@Override
	public boolean isConcised() {
		return false;
	}

	@Override
	public boolean isStretchable() {
		return isVisible();
	}
	
	@Override
    public Integer getHeight() {
        if (isVisible()) {
            return 12;
        } else {
            return 0;
        }
    }
	
	public void changeSNPView(Class<? extends AreaRequestHandler> handler) {
        // turn on highlighting mode
        changeView = true;
        
        // set reference data
        view.getQueueManager().createQueue(file, handler);
        view.getQueueManager().addResultListener(file, this);
    }
    
    public void returnSNPView() {
        // turn off highlighting mode
        changeView = false;
    }
    
    @Override
    public String getName() {
    	return "SNPTrack";
    }
	
}
