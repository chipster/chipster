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
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SNPTrack extends Track {
	Collection<RegionContent> values = new LinkedList<RegionContent>();
	Color a = new Color(64, 192, 64, 128);
	Color c = new Color(64, 64, 192, 128);
	Color g = new Color(128, 128, 128, 128);
	Color t = new Color(192, 64, 64, 128);
	Color frameShift = new Color(255, 255, 255, 255);

	boolean changeView = false;
	
	enum ConsequenceName {
		NONE, UPSTREAM
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
	            String allele = (String)value.values.get(ColumnType.ALLELE);
	            long position = getView().bpToTrack(value.region.start);
	            if (changeView) {

	            } else {
	            	//deal with frameshift
	            	if (allele.startsWith("A")) {
		            	drawables.add(new RectDrawable((int)position, 1, 3, 5, a, a));
		            } else if (allele.startsWith("C")) {
		            	drawables.add(new RectDrawable((int)position, 1, 3, 5, c, c));
		            } else if (allele.startsWith("G")) {
		            	drawables.add(new RectDrawable((int)position, 1, 3, 5, g, g));
		            } else if (allele.startsWith("T")) {
		            	drawables.add(new RectDrawable((int)position, 1, 3, 5, t, t));
		            }
	            }
	            
			}
		}
		
		return drawables;
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
            return 5;
        } else {
            return 0;
        }
    }
	
}
