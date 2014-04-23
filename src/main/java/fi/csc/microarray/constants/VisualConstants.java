package fi.csc.microarray.constants;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Stroke;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.plaf.BorderUIResource;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;
import com.jgoodies.looks.plastic.Plastic3DLookAndFeel;
import com.jgoodies.looks.plastic.theme.ExperienceBlue;

import fi.csc.microarray.client.Session;

public class VisualConstants {
	
    // icons from Eclipse
    public static final ImageIcon PHENODATA_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/phenodata.gif"));
    public static final ImageIcon IMPORT_CANCEL_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/import_cancel.gif"));
    public static final ImageIcon IMPORT_NEXT_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/import_forward.gif"));
    public static final ImageIcon IMPORT_FINISH_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/import_finish.png"));
    public static final ImageIcon IMPORT_BACK_ICON = 
       new ImageIcon(VisualConstants.class.getResource("/eclipse/import_backward.gif"));
    public static final ImageIcon IMPORT_RESET_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/eclipse/import_reset.png"));
    public static final ImageIcon IMPORT_HEADER_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/eclipse/import_header.png"));
    public static final ImageIcon IMPORT_FOOTER_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/eclipse/import_footer.png"));
    public static final ImageIcon IMPORT_TITLE_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/eclipse/import_title.png")); 
    public static final ImageIcon DEFAULT_VIEW_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/defaultView.gif"));
    public static final ImageIcon EXPORT_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/export.gif"));
    public static final ImageIcon PASTE_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/paste.gif"));
    public static final ImageIcon COPY_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/copy.gif"));
    public static final ImageIcon DELETE_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/delete.gif"));
    public static final ImageIcon HELP_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/help.gif"));
    public static final ImageIcon LINK_PHENODATA_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/link_phenodata.png"));
    public static final ImageIcon LINK_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/link.gif"));
    public static final ImageIcon UNLINK_MENUICON = 
        new ImageIcon(VisualConstants.class.getResource("/eclipse/unlink.gif"));
    public static final ImageIcon QUICKLINK_ICON =
        new ImageIcon(VisualConstants.class.getResource("/eclipse/bulb.gif"));

	// rest of the icons done by us (mostly by Mikko Koski and Petri Klemelä)
	public static ImageIcon ZOOM_IN_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/zoom-in.png"));
	public static ImageIcon ZOOM_OUT_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/zoom-out.png"));
	public static ImageIcon MAGNIFIER_ICON = 
	    new ImageIcon(VisualConstants.class.getResource("/viewmag.png"));
    public static ImageIcon CLOSE_FILE_ICON = 
	    new ImageIcon(VisualConstants.class.getResource("/fileclose.png"));
	public static final ImageIcon APPLICATION_ICON =
        new ImageIcon(VisualConstants.class.getResource("/chipster_icon.png"));
	public static final ImageIcon LOGIN_BANNER =
        new ImageIcon(VisualConstants.class.getResource("/login_banner.png"));
    public static final ImageIcon QUERY_PARAM_CLASS_ICON =
        new ImageIcon(VisualConstants.class.getResource("/circle.gif"));
    public static final ImageIcon QUERY_PARAM_INSTANCE_ICON =
        new ImageIcon(VisualConstants.class.getResource("/square.gif"));    

    public static ImageIcon ZOOM_IN_CURSOR_IMAGE= 
		new ImageIcon(VisualConstants.class.getResource("/zoomInCursor.png"));
	public static ImageIcon ZOOM_OUT_CURSOR_IMAGE = 
		new ImageIcon(VisualConstants.class.getResource("/zoomOutCursor.png"));
	public static ImageIcon ROTATE_CURSOR_IMAGE = 
		new ImageIcon(VisualConstants.class.getResource("/3dRotateCursor.png"));
	public static ImageIcon ROTATE_AND_ZOOM_CURSOR_IMAGE = 
		new ImageIcon(VisualConstants.class.getResource("/rotatezoom.png"));
	public static ImageIcon ROTATE_IMAGE = 
		new ImageIcon(VisualConstants.class.getResource("/3dRotate.png"));
	public static ImageIcon ARROW_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/arrow.png"));
	public static ImageIcon HAND_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/hand.png"));
	public static ImageIcon SHOW_ALL_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/showall.png"));
	public static final Icon DOUBLE_FORWARD_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/forward.png"));		
	public static final Icon DOUBLE_FORWARD_BW_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/forward_bw.png"));
	public static final ImageIcon SUITABLE_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/yes.png"));
    public static final ImageIcon INCOMPATIBLE_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/no.png"));    
    public static final ImageIcon SUITABILITY_WARNING_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/maybe.png"));
    public static ImageIcon GENERATE_HISTORY_ICON =    	
		new ImageIcon(VisualConstants.class.getResource("/history.png"));
    public static ImageIcon TO_TOP_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/toTop.png"));
    public static ImageIcon TO_BOTTOM_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/toBottom.png"));
    public static ImageIcon REDRAW_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/redraw_bw.png"));    
    public static ImageIcon MAXIMISE_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/maximise.png"));
    public static ImageIcon RESTORE_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/restore.png"));
    public static ImageIcon QUESTION_MARK_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/question_mark.png"));
    public static ImageIcon SPLIT_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/split.png"));
    public static ImageIcon TO_WINDOW_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/to_window.png"));
    public static ImageIcon CLOSE_ICON = 
		new ImageIcon(VisualConstants.class.getResource("/close.png"));
    
    public static ImageIcon XY_PLANE = 
		new ImageIcon(VisualConstants.class.getResource("/XYPlane.png"));
    public static ImageIcon XZ_PLANE = 
		new ImageIcon(VisualConstants.class.getResource("/XZPlane.png"));
    public static ImageIcon YZ_PLANE = 
		new ImageIcon(VisualConstants.class.getResource("/YZPlane.png"));

    
//Crystal icons, LGPL, http://www.everaldo.com/crystal/
    public static ImageIcon PHENODATA_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/kate.png"));
    public static ImageIcon TEXT_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/txt.png"));
    public static ImageIcon IMAGE_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/no3D.png"));
    public static ImageIcon HTML_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/agt_web.png"));
    
//Modified icons based on the Crystal icons, LGPL, http://www.everaldo.com/crystal/   
    public static ImageIcon PDF_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/modified/man.png"));
    public static ImageIcon EXT_BROWSER_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/modified/ext-browser.png"));
    public static ImageIcon SPREADSHEET_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/modified/spreadsheet.png"));
    public static ImageIcon ARRAY_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/crystal-LGPL/modified/array.png"));
    
// rest of the icons done by us (mostly by Mikko Koski and Petri Klemelä)
    public static ImageIcon VENN_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/venndi-48.png"));
    public static ImageIcon EMPTY_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/empty.png"));
    public static ImageIcon SOM_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/som2.png"));
    public static ImageIcon SCATTER_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/scatter4.png"));
    public static ImageIcon PROFILE_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/profile2.png"));
    public static ImageIcon SCATTER3DPCA_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/pca.png"));
    public static ImageIcon HC_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/hc3.png"));
    public static ImageIcon SCATTER3D_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/scatter3d_2.png"));
    public static ImageIcon PROFILES_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/profiles3.png"));
    public static ImageIcon VOLCANO_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/volcano2.png"));
    public static ImageIcon HISTOGRAM_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/histogram2.png"));
    public static ImageIcon GB_MENUICON = 
    		new ImageIcon(VisualConstants.class.getResource("/listicons/gb.png"));
    
    
    public static final ImageIcon SPLASH_SCREEN =
        new ImageIcon(VisualConstants.class.getResource("/splash.png"));
    public static final Color SPLASH_BORDER_COLOR =
    	new Color(150, 150, 150);
    public static final Font SPLASH_SCREEN_FONT = 
    	new Font("SansSerif", Font.PLAIN, 10);
	
    public static final Font VISUALISATION_TITLE_FONT = 
    	new Font("SansSerif", Font.PLAIN, 16);
    
    public static final Font MONOSPACED_FONT = 
    	new Font("Monospaced", Font.PLAIN, 11);



    public static final ImageIcon ICON_TYPE_BINARY =
        new ImageIcon(VisualConstants.class.getResource("/types/binary.png"));
    public static final ImageIcon ICON_TYPE_FOLDER =
        new ImageIcon(VisualConstants.class.getResource("/types/folder.gif"));
    public static final ImageIcon ICON_TYPE_HTML =
        new ImageIcon(VisualConstants.class.getResource("/types/html.png"));
    public static final ImageIcon ICON_TYPE_IMAGE =
        new ImageIcon(VisualConstants.class.getResource("/types/image.png"));
    public static final ImageIcon ICON_TYPE_PHENODATA =
        new ImageIcon(VisualConstants.class.getResource("/types/phenodata.png"));
    public static final ImageIcon ICON_TYPE_RAWDATA =
        new ImageIcon(VisualConstants.class.getResource("/types/rawdata.png"));
    public static final ImageIcon ICON_TYPE_TABLE =
        new ImageIcon(VisualConstants.class.getResource("/types/table.png"));
    public static final ImageIcon ICON_TYPE_TEXT =
        new ImageIcon(VisualConstants.class.getResource("/types/text.png"));

    public static final ImageIcon ARROW_UP_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/arrow_up.png"));
    public static final ImageIcon ARROW_DOWN_ICON = 
    	new ImageIcon(VisualConstants.class.getResource("/arrow_down.png"));
    public static final ImageIcon STOP_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/stop.png"));    
    public static final ImageIcon RUNNING_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/running.gif"));
    
    public static final ImageIcon OPEN_SESSION_LINK_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/tree.png"));
    public static final ImageIcon IMPORT_LINK_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/table.png"));
    public static final ImageIcon EXAMPLE_SESSION_ICON = 
        new ImageIcon(VisualConstants.class.getResource("/try.png"));
    

    public static final int LEFT_PANEL_WIDTH = 400;
    public static final int TREE_PANEL_HEIGHT = 270;
    public static final int GRAPH_PANEL_HEIGHT = 240;
        
    // for some child screens
    public static final Dimension DEFAULT_SCREEN_DIMENSION = new Dimension(640,480); 
    
    // colors for Plastic3D look and feel
    public static final Color PLASTIC3D_FOCUS_COLOR = new Color(200, 200, 200);
    public static final Color TOOL_LIST_BORDER_COLOR = new Color(128,128,128);
    public static final Color PHENODATA_TABLE_UNEDITABLE_CELL_BACKGROUND = new Color(245, 245, 245);
    
    public static final float DEFAULT_FONT_SIZE = 12f;
    /**
     * Color for uneditable JTextArea background. Used in details pane and 
     * ErrorScreen for example. 
     */
	public static Color TEXTAREA_UNEDITABLE_BACKGROUND = UIManager.getColor("Panel.background");
    
    public static UIDefaults getUIDefaults() {
    	
    	//Uncomment to see possible key values
    	//listUIDefaults();
    	
    	UIDefaults defaults = new UIDefaults();

    	// The defauls are specified for each look n feel and color theme
    	if(UIManager.getLookAndFeel() instanceof Plastic3DLookAndFeel 
    			&& Plastic3DLookAndFeel.getPlasticTheme() instanceof ExperienceBlue){
    		
    		BorderUIResource emptyBorder = new BorderUIResource(BorderFactory.createEmptyBorder());
    		
    		// Removes borders from menubar
    		defaults.put("MenuBar.border", emptyBorder);
    		
    		defaults.put("SplitPaneDivider.border", emptyBorder);
    		defaults.put("SplitPane.border", emptyBorder);
    		defaults.put(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);
    		
    		defaults.put("TaskPane.titleBackgroundGradientStart", 
    				UIManager.getColor("Panel.background"));
    		defaults.put("TaskPane.titleBackgroundGradientEnd", 
    				UIManager.getColor("Panel.background"));
    		defaults.put("TaskPaneContainer.background", Color.white);
    		
    		defaults.put("SimpleInternalFrame.activeTitleForeground", Color.white);
    		
    		// Adds textarea background to white. This affects for example 
    		// help textarea on the top right corner and affymetrix wizard
    		defaults.put("TextArea.background", Color.WHITE);    		
    		
    		TEXTAREA_UNEDITABLE_BACKGROUND = UIManager.getColor("Panel.background");
    	} else {
    		// There is no specified look and feel options for this LAF. Use defaults.
    	}
    	
    	return defaults;
    }
    
    /**
     * Prints list of UIDefault keys which are currently in use. Note that 
     * this is not the whole list of available keys. As we know, there is no 
     * documented list of all keys. Keys that are not in this list could 
     * be found by looking source codes of Look And Feel classes, 
     * javax.swing.plaf.metal.MetalLookAndFeel for example.
     * 
     * @see javax.swing.plaf.metal.MetalLookAndFeel
     *
     */
    public static void listUIDefaults(){
    	ArrayList<String> defaultKeys = new ArrayList<String>();
    	for (Enumeration<Object> e = UIManager.getDefaults().keys(); e.hasMoreElements(); ){
    		Object obj = e.nextElement();
    		if (obj instanceof String){
    			defaultKeys.add(obj.toString());
    		}
    	}
    	String[] keysAsArray = defaultKeys.toArray(new String[defaultKeys.size()]);
    	Arrays.sort(keysAsArray);
    	
    	for(String s: keysAsArray){
    		System.out.println(s);
    	}
    }
    
    /**
     * See docs/category_color_scheme.jpg for color scheme.
     */
    public final static Color[] CATEGORY_COLORS = new Color[]{
		new Color(195,182,162), //used for raw data
		new Color(213,199,150), 
		new Color(231, 223, 112), 
		new Color(213,159,69),
		new Color(231,136,28),
		new Color(213,56,51),
		new Color(128,163,183),
		new Color(1,119,183),
		new Color(98, 154, 155),
		new Color(164,153,0),
		new Color(131,1,11),
		new Color(192,210,222),
		new Color(71,139,140)
	};
	
	//Used for the selection rectangles in Scatterplots
	public static final Stroke dashLine = 
		new BasicStroke(
				1, 
				BasicStroke.CAP_BUTT,
				BasicStroke.JOIN_BEVEL,
				0,
				new float[] {9}, 
				0
		);

	public static final String HTML_DIALOG_TITLE_STYLE = "\"font-weight:bold;font-size:115%\"";
	
	public static final String getImportDirectlyText() {
		return Session.getSession().getApplication().isStandalone() ? "Import directly" : "Import directly if possible"; 
	}

	public static final Color COLOR_BEIGE = new Color(0xc3b6a2);
	public static final Color COLOR_SAND = new Color(0xd5c796);
	public static final Color COLOR_YELLOW = new Color(0xe7df70);
	public static final Color COLOR_BROWN_LIGHT = new Color(0xd59f45);
	public static final Color COLOR_ORANGE = new Color(0xe7881c);
	public static final Color COLOR_RED = new Color(0xd53833);
	public static final Color COLOR_BLUE_GREY = new Color(0x80a3b7);
	public static final Color COLOR_BLUE = new Color(0x0177b7);
	public static final Color COLOR_BLUE_EVEN_BRIGHTER = new Color(0x01A6FF);
	public static final Color COLOR_BLUE_BRIGHTER = new Color(0x0199EB);
	public static final Color COLOR_BLUE_GREEN = new Color(0x629a9b);
	public static final Color COLOR_YELLOW_DIRTY = new Color(0xa49900);
	public static final Color COLOR_BROWN = new Color(0x83010b);
	public static final Color COLOR_BLUE_LIGHT = new Color(0xc0d2de);
}
