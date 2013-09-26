package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;

public class RepeatMaskerTrack extends Track{

	private Color color;

	private Collection<Feature> regions = new TreeSet<Feature>();

	public RepeatMaskerTrack(){
		super();
		this.color = Color.lightGray;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();      

		if (regions != null) {

			Iterator<Feature> iter = regions.iterator();

			Feature regionContent = null;

			while (iter.hasNext()) {

				regionContent = iter.next();

				if (!getView().requestIntersects(regionContent.region)) {
					iter.remove();
					continue;
				}

				BpCoord startBp = regionContent.region.start;
				//Increase by one to fill the last base pair
				BpCoord endBp = new BpCoord(regionContent.region.end.bp + 1, regionContent.region.end.chr);
				

				long startX = getView().bpToTrack(startBp);
				long endX = getView().bpToTrack(endBp);

				drawables.add(new RectDrawable(new Rectangle((int)startX, 0,
						(int)(endX-startX), this.getTrackHeight()), color, color));
			}
		}


		return drawables;
	}

	@Override
	public void processDataResult(DataResult dataResult) {

		this.regions.addAll(dataResult.getFeatures());
	}

	@Override
	public int getTrackHeight() {
		return 5;
	}
	
    @Override
	public void defineDataTypes() {
    	addDataType(DataType.REGION);
	}
	
    @Override
    public String getTrackName() {
    	return "RepeatMaskerTrack";
    }
}
