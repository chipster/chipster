package fi.csc.microarray.client.visualisation.methods;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.entity.ChartEntity;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.PlotState;
import org.jfree.util.ObjectUtilities;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.databeans.DataBean;

public class VenndiPlot extends Plot implements ChartMouseListener, Cloneable, Serializable {
	
	public class VennDataset {
		private String[][] datas;
		private List<DataBean> datasets;
		private Map<DataBean, Map<String, Integer>> indexMaps;
		
		
		/**
		 * Two dimensional table of indetifiers. First dimension has to be 7 and obey the order
		 * listed in enum VenndiPlot.AREAS. Other dimension lists all the identifiers for those
		 * areas.
		 * 
		 * @param data
		 * @param ids 
		 */
		public VennDataset (String[][] data, List<DataBean> datas, Map<DataBean, Map<String, Integer>> indexMaps2){
			this.datas = data;
			this.datasets = datas;
			this.indexMaps = indexMaps2;
		}

		public int getCount(VenndiPlot.AREAS area){
			return datas[area.ordinal()].length; 
		}
		
		public String getDatasetName(AREAS area) {
			if( area.ordinal() < datasets.size()){
				return datasets.get(area.ordinal()).getName();
			} 
			return "";
		}
		
		public List<DataBean> getDataBeans(){
			return datasets;
		}

		public Collection<String> getIdentifiers(AREAS areas) {
			return Arrays.asList(datas[areas.ordinal()]);
		}
		
		public Map<DataBean, Set<Integer>> getIndexes(AREAS area){
			
			Map<DataBean, Set<Integer>> indexes = new HashMap<DataBean, Set<Integer>>();
			
			for(DataBean data : datasets){
				indexes.put(data, new HashSet<Integer>());
				
				for(String id : getIdentifiers(area)){					
					indexes.get(data).add(indexMaps.get(data).get(id));
				}
			}
			
			return indexes;
		}
	}
    
	private static final long serialVersionUID = 1L;

	/** The cell information text is drawn with this font. */
    private Font descriptionFont;

	private VennDataset dataset;

	private Area[] areas = new Area[7];
	private Set<Area> selected = new HashSet<Area>();
	
	private VennDiagram visualisation;
     
    /** The resourceBundle for the localization. */
    protected static ResourceBundle localizationResources 
        = ResourceBundle.getBundle("org.jfree.chart.plot.LocalizationBundle");
    
    
    public VenndiPlot(String[][] idTable, List<DataBean> datas, Map<DataBean, Map<String, Integer>> indexMaps, VennDiagram parent) throws NullPointerException {
    	this.dataset = new VennDataset(idTable, datas, indexMaps);
    	
        this.visualisation = parent;
        
        //Fill variable areas to be able to mark selections from the application
        updateAreas(new Rectangle(), null);
    }
    
    public List<String> getSelected() {
    	
    	List<String> selectedIds = new ArrayList<String>();
    	
    	for (int i = 0; i < areas.length; i++){
    		if (selected.contains(areas[i])){
    			selectedIds.addAll(dataset.getIdentifiers(AREAS.values()[i]));
    		}
    	}
    	
    	return selectedIds;
    }
    
    public void setSelected(Collection<String> selectedIds){
    	selected.clear();
    	
    	if(selectedIds != null){
    		for (AREAS area : AREAS.values()){
    			if (selectedIds.containsAll(dataset.getIdentifiers(area)) && 
    					dataset.getCount(area) > 0){
    				selected.add(areas[area.ordinal()]);
    			}
    		}
    	}
    }
    
    
    public VennDataset getDataset() {
        return this.dataset;
    }

        
    public String getPlotType() {
        //return localizationResources.getString("VennDiagram");
    	return "VennDiagram";
    }

    /**
     * Returns the font used for cell descriptions.
     *
     * @return Description font.
     */
    public Font getDescriptionFont() {
        return this.descriptionFont;
    }

    public enum AREAS { A, B, C, AB, AC, BC, ABC };
        
    /**
     * Draws the plot on a Java2D graphics device (such as the screen or 
     * a printer).
     *
     * @param g2  the graphics device.
     * @param area  the area within which the plot should be drawn.
     * @param anchor  the anchor point (<code>null</code> permitted).
     * @param parentState  the state from the parent plot, if there is one
     * (<code>null</code> permitted.)
     * @param info  collects info about the drawing (<code>null</code> permitted).
     * @throws NullPointerException  if g2 or area is null.
     */ 
    @Override
    public void draw(java.awt.Graphics2D g2, java.awt.geom.Rectangle2D area, 
            java.awt.geom.Point2D anchor, PlotState parentState, 
            PlotRenderingInfo info) {
    	    	
    	
        // adjust for insets...
        this.getInsets().trim(area);

        if (info != null) {
            info.setPlotArea(area);
            info.setDataArea(area);
        }

        //this.setBackgroundPaint(Color.black);
        
        drawBackground(g2, area);
        drawOutline(g2, area);

        Shape savedClip = g2.getClip();
        g2.clip(area);
/*
        Composite originalComposite = g2.getComposite();
        g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 
                getForegroundAlpha()));
*/
        drawVenn(g2, area, info);
        
        g2.setClip(savedClip);
   //   g2.setComposite(originalComposite);

        drawOutline(g2, area);
    }          
    
    private void drawVenn(Graphics2D g2, java.awt.geom.Rectangle2D area, 
            PlotRenderingInfo info) {            	
    		
    	updateAreas(area, info);
    	
   	      	Color[] colors = {
   			VisualConstants.CATEGORY_COLORS[4],
   			VisualConstants.CATEGORY_COLORS[5],
   			VisualConstants.CATEGORY_COLORS[7],
   			new Color(222, 98, 39),
   			new Color(115, 127, 107),
   			new Color(107, 88, 117),
   			VisualConstants.CATEGORY_COLORS[10] 
   	};

    	
    	
//    	  //Little bit different colors
//    	  	Color[] colors = {
//    			VisualConstants.CATEGORY_COLORS[4],
//    			VisualConstants.CATEGORY_COLORS[5],
//    			VisualConstants.CATEGORY_COLORS[7],
//    			VisualConstants.CATEGORY_COLORS[3],
//    			VisualConstants.CATEGORY_COLORS[6],
//    			new Color(147, 147, 147),
//    			VisualConstants.CATEGORY_COLORS[10] 
//    	};
    	  
    	    	
   	     g2.setFont(g2.getFont().deriveFont(g2.getFont().getSize2D() * 0.8f));
    	
    	for ( int i = 0; i < areas.length; i++){
    		Shape shape = areas[i];
    		if(shape != null){
    			if(selected.contains(shape)){
    				
    				Color c = new Color(colors[i].getRed(), colors[i].getGreen(), colors[i].getBlue(), 200);
    	
    				g2.setPaint(c);
    				g2.fill(shape);
    				g2.setPaint(Color.black);
    				g2.setStroke(new BasicStroke(3));
    				g2.draw(shape);
    				g2.setFont(g2.getFont().deriveFont(Font.BOLD));

    			} else {
    				
    				Color c = new Color(colors[i].getRed(), colors[i].getGreen(), colors[i].getBlue(), 255);
    	
    				g2.setPaint(c);
    				g2.fill(shape);    				
    				g2.setPaint(Color.black);
    				g2.setStroke(new BasicStroke(1));
    				g2.draw(shape);
    				g2.setFont(g2.getFont().deriveFont(Font.PLAIN));
    			}
    			
    			g2.setPaint(Color.black);
    			String str = "" + dataset.getCount(AREAS.values()[i]);
    			
    			int x = (int)shape.getBounds().getCenterX() - 
    				SwingUtilities.computeStringWidth(g2.getFontMetrics(), str) / 2;
    			
    			int y = (int)shape.getBounds().getCenterY() + g2.getFontMetrics().getHeight() / 2;
    			
    			g2.drawString(str, x, y);
    			
    			str = dataset.getDatasetName(AREAS.values()[i]);
    			
    			x = (int)shape.getBounds().getCenterX() - 
				SwingUtilities.computeStringWidth(g2.getFontMetrics(), str) / 2;
    			
    			y += g2.getFontMetrics().getHeight();
    			
    			g2.drawString(str, x, y);
    		}    		    		
    	}    
    }
    
    private void updateAreas(Rectangle2D area, PlotRenderingInfo info){

    	double r = area.getHeight() < area.getWidth() ? (int)area.getHeight()/4 : (int)area.getWidth()/4;
    	
    	//Distance of circle center from the area center		
    	double d = r * 0.7;
    	
    	//Components of the d in angle of 30 degrees
    	double dx = d * Math.cos(Math.PI / 6);
    	double dy = d * Math.sin(Math.PI / 6);    	    	
    	
    	Area[] circles = new Area[3];		
    	circles[AREAS.A.ordinal()] = new Area(new Ellipse2D.Double(area.getCenterX() - r, 		area.getCenterY() - d - r, 2*r, 2*r));
    	circles[AREAS.B.ordinal()] = new Area(new Ellipse2D.Double(area.getCenterX() - dx - r, 	area.getCenterY() + dy - r, 2*r, 2*r));
    	circles[AREAS.C.ordinal()] = new Area(new Ellipse2D.Double(area.getCenterX() + dx - r, 	area.getCenterY() + dy - r, 2*r, 2*r));
    	
    	if(areas[0] == null){
    		for (int i = 0; i < areas.length; i++){
    			areas[i] = new Area();
    		}
    	} else {
    		for (Area a : areas){
    			a.reset();
    		}
    	}

    	areas[AREAS.A.ordinal()].add(circles[AREAS.A.ordinal()]);
    	areas[AREAS.A.ordinal()].subtract(circles[AREAS.B.ordinal()]);
    	areas[AREAS.A.ordinal()].subtract(circles[AREAS.C.ordinal()]);
    	
    	areas[AREAS.B.ordinal()].add(circles[AREAS.B.ordinal()]);
    	areas[AREAS.B.ordinal()].subtract(circles[AREAS.A.ordinal()]);
    	areas[AREAS.B.ordinal()].subtract(circles[AREAS.C.ordinal()]);
    	
    	areas[AREAS.C.ordinal()].add(circles[AREAS.C.ordinal()]);
    	areas[AREAS.C.ordinal()].subtract(circles[AREAS.A.ordinal()]);
    	areas[AREAS.C.ordinal()].subtract(circles[AREAS.B.ordinal()]);
    	
    	areas[AREAS.AB.ordinal()].add(circles[AREAS.A.ordinal()]);
    	areas[AREAS.AB.ordinal()].intersect(circles[AREAS.B.ordinal()]);
    	areas[AREAS.AB.ordinal()].subtract(circles[AREAS.C.ordinal()]);
    	
    	areas[AREAS.BC.ordinal()].add(circles[AREAS.B.ordinal()]);
    	areas[AREAS.BC.ordinal()].subtract(circles[AREAS.A.ordinal()]);
    	areas[AREAS.BC.ordinal()].intersect(circles[AREAS.C.ordinal()]);
    	
    	areas[AREAS.AC.ordinal()].add(circles[AREAS.C.ordinal()]);
    	areas[AREAS.AC.ordinal()].intersect(circles[AREAS.A.ordinal()]);
    	areas[AREAS.AC.ordinal()].subtract(circles[AREAS.B.ordinal()]);
    	
    	areas[AREAS.ABC.ordinal()].add(circles[AREAS.A.ordinal()]);
    	areas[AREAS.ABC.ordinal()].intersect(circles[AREAS.B.ordinal()]);
    	areas[AREAS.ABC.ordinal()].intersect(circles[AREAS.C.ordinal()]);
    	
    	if(info != null){
    		for(Area shape : areas){
    			createEntity(shape, info);
    		}
    	}
    }
    
    /**
     * Implements the ChartMouseListener interface. This
     * method does nothing.
     *
     * @param event  the mouse event.
     */ 
    public void chartMouseMoved(ChartMouseEvent event) {
        ;
    }
    
    public void chartMouseClicked(ChartMouseEvent e) {    	    	                
    	
    	if(e.getEntity() != null){
    		for(Area shape : areas){    		
    			//if(shape.contains(e.getTrigger().getPoint())){
    			if(shape.equals(e.getEntity().getArea())){

    				if(!e.getTrigger().isControlDown()){
    					selected.clear();
    				}
    				if(selected.contains(shape)){
    					selected.remove(shape);
    				} else {
    					selected.add(shape);
    				}
    			}    	
    		}
    	} else {    		
    		selected.clear();
    	}
    	    	
    	setSelectedForList(true);
    }
    
	public void setSelectedForList(boolean dispatchEvent){
		
    	List<String> ids = new ArrayList<String>();
    	// The rowSelectionManager needs row indexes
    	Map<DataBean, Set<Integer>> indexes = new HashMap<DataBean, Set<Integer>>();
    	
    	// Due to multiple datas the row indexes are meaningful only within one specific data
    	for(int i = 0; i < areas.length; i++){
    		if(selected.contains(areas[i])){    		
    			ids.addAll(dataset.getIdentifiers(AREAS.values()[i]));
    			Map<DataBean, Set<Integer>> map = dataset.getIndexes(AREAS.values()[i]);
    			for(DataBean data: map.keySet()){
    				if(!indexes.containsKey(data)){
    					indexes.put(data, new HashSet<Integer>());
    				}
    				
    				indexes.get(data).addAll(map.get(data));
    			}
    		}
    	}
		
		visualisation.setSelectedListContent(ids, indexes, this, dispatchEvent);
	}
    
    /**
     * Tests this plot for equality with an arbitrary object.  Note that the 
     * plot's dataset is NOT included in the test for equality.
     *
     * @param obj  the object to test against (<code>null</code> permitted).
     *
     * @return <code>true</code> or <code>false</code>.
     */
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof VenndiPlot)) {
            return false;
        }
        VenndiPlot that = (VenndiPlot) obj;

//        if (!ObjectUtilities.equal(this.urlGenerator, that.urlGenerator)) {
//            return false;
//        }
        
        if (!ObjectUtilities.equal(this.descriptionFont, 
                that.descriptionFont)) {
            return false;
        }
    
        // can't find any difference...
        return true;
    }

    /**
     * Provides serialization support.
     *
     * @param stream  the output stream.
     *
     * @throws IOException  if there is an I/O error.
     * @throws NullPointerException  if stream is null.
     */
    private void writeObject(ObjectOutputStream stream) throws IOException {
        stream.defaultWriteObject();
    }

    /**
     * Provides serialization support.
     *
     * @param stream  the input stream.
     *
     * @throws IOException  if there is an I/O error.
     * @throws ClassNotFoundException  if there is a classpath problem.
     * @throws NullPointerException  if stream is null.
     */
    private void readObject(ObjectInputStream stream) 
        throws IOException, ClassNotFoundException {
        stream.defaultReadObject();
    }
    
    private void createEntity(Shape shape, PlotRenderingInfo info) {
        
        if (info != null) {
            EntityCollection entities = info.getOwner().getEntityCollection();
            if (entities != null) {
                
                entities.add(
                        new ChartEntity(shape));
            }
        }
    }
}

