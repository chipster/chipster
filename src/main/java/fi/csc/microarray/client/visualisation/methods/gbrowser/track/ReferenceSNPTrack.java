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
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * <p>Track showing known SNP's of the genome.</p> 
 * 
 * <p>There are two ways of showing data: 1. by nucleotide color, or 2. by consequence to transcript.</p>
 * 
 * <p>In the consequence to transcript way the color is get by simple prioritization
 * mechanism. There is array of colors, which are ordered by importance and
 * consequences enumeration, which are also ordered by importance, and particular
 * color refer to particular consequence name. If consequence name number is greater than
 * maxColorIndex, then 'default' gray color is defined.</p>
 * 
 * <p>If importance sequence is changed, so the colors array must be changed accordingly.</p>
 * 
 * @author Vilius Zukauskas
 *
 */

public class ReferenceSNPTrack extends Track {
	
	private static int width = 8;
	Collection<RegionContent> values = new LinkedList<RegionContent>();
	Long lastPosition;
	long maxBpLength;
	long minBpLength;
	Color a = new Color(64, 192, 64, 128); //green
	Color c = new Color(64, 64, 192, 128); //blue
	Color g = new Color(128, 128, 128, 128); // gray
	Color t = new Color(192, 64, 64, 128); //red
	Color forExceptions = new Color(139, 69, 19, 128);//brown
	
	//ordered by prioritiazation
	private Color[] colors = new Color[] {
			new Color(255, 0, 0, 128),//red
			new Color(255, 0, 0, 128),//red
			new Color(255, 105, 180, 128),//pink
			new Color(255, 215, 0, 128),//gold
			new Color(50, 205, 50, 128),//lemon green
			new Color(0, 0, 255, 128), //blue
			new Color(190, 190, 190, 128),//gray
	};
	
	private int maxColorIndex = colors.length-1;
	
	private int lastConsequenceNumber = -1;
	
	boolean changeView = true;
	
	//ordered according colors
	enum ConsequenceName {
		STOP_GAINED, STOP_LOST, FRAMESHIFT_CODING, SYNONYMOUS_CODING,
		NON_SYNONYMOUS_CODING, WITHIN_MATURE_miRNA, UPSTREAM, DOWNSTREAM, INTRONIC, INTERGENETIC,  
		COMPLEX_INDEL, PARTIAL_CODON, REGULATORY_REGION, WITHIN_MATURE_mIRNA,
		PRIME5_UTR,	PRIME3_UTR , STOP_GAINED_FRAMESHIFT, STOP_GAINED_SPLICE_SITE,
		STOP_LOST_SPLICE_SITE, FRAMESHIFT_CODING_SPLICE_SITE, 
		STOP_GAINED_FRAMESHIFT_CODING_SPLICE_SITE, NON_SYNONYMOUS_CODING_SPLICE_SITE,
		SPLICE_SITE_SYNONYMOUS_CODING, SPLICE_SITE_5PRIME_UTR, SPLICE_SITE_3PRIME_UTR, 
		ESSENTIAL_SPLICE_SITE_INTRONIC, SPLICE_SITE_INTRONIC, WITHIN_NON_CODING_GENE, 
		INTRONIC_NMD_TRANSCRIPT, NONE,
	}

	public ReferenceSNPTrack(View view, DataSource file, long minBpLength, long maxBpLength) {
		super(view, file);
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		Iterator<RegionContent> iter = values.iterator();
		
		if (values != null) {
			while (iter.hasNext()) {
				RegionContent value = iter.next();
				
				// remove those that are not in this view
	            if (!value.region.intersects(getView().getBpRegion())) {
	                iter.remove();
	                continue;
	            }
	            
	            if (lastPosition == null) {
	            	lastPosition = (long)getView().bpToTrack(value.region.start);
	            }
	            
	            if (lastConsequenceNumber == -1) {
					lastConsequenceNumber = ConsequenceName.valueOf(
							(String)value.values.get(ColumnType.CONSEQUENCE_TO_TRANSCRIPT)).ordinal();
				}
	            String allele = (String)value.values.get(ColumnType.ALLELE);
	            long position = getView().bpToTrack(value.region.start);
	            if (changeView) {
	            	String consequence = (String)value.values.get(ColumnType.CONSEQUENCE_TO_TRANSCRIPT);
	            	
        			int now;
					try {
						now = ConsequenceName.valueOf(consequence).ordinal();
					} catch (Exception e) {
						//if it's a consequence we don't know about
						now = ConsequenceName.NONE.ordinal();
					}
        			if (lastPosition != position) {
        				drawables.add(new RectDrawable((int)position, 1, 
	            				width, getHeight(), colors[(now > maxColorIndex) ? maxColorIndex : now],
	            				colors[(now > maxColorIndex) ? maxColorIndex : now]));
        			}
	            	if (lastConsequenceNumber < now) {
	            		drawables.add(new RectDrawable((int)position, 1, 
	            				width, getHeight(), colors[(now > maxColorIndex) ? maxColorIndex : now],
	            				colors[(now > maxColorIndex) ? maxColorIndex : now]));
	            	} else {
	            		drawables.add(new RectDrawable((int)position, 1, 
	            				width, getHeight(), colors[(lastConsequenceNumber > maxColorIndex) ? maxColorIndex : lastConsequenceNumber],
	            				colors[(lastConsequenceNumber > maxColorIndex) ? maxColorIndex : lastConsequenceNumber]));
	            	}

	            } else {
            		if (allele.matches("[ACGT]/A")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), a, a));
		            } else if (allele.matches("[ACGT]/C")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), c, c));
		            } else if (allele.matches("[ACGT]/G")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), g, g));
		            } else if (allele.matches("[ACGT]/T")) {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), t, t));
		            } else {
		            	drawables.add(new RectDrawable((int)position, 1, width, getHeight(), forExceptions, forExceptions));
	            	}	            		
	            }
	            lastPosition = position;
	            lastConsequenceNumber = ConsequenceName.valueOf(
						(String)value.values.get(ColumnType.CONSEQUENCE_TO_TRANSCRIPT)).ordinal();
			}
		}

		return drawables;
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {
		for (RegionContent content : areaResult.getContents()) {

			if (content.values.get(ColumnType.STRAND) == getStrand()) {
				values.add(content);
			}
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

	public void changeSNPView() {
        // turn on highlighting mode
        changeView = true;
        
        // set reference data
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
    
    @Override
    public boolean isVisible() {
    	if (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength) {
    		return true;
    	} else {
    		return false;
    	}
    }

}
