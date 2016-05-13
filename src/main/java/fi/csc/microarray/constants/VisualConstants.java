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
	public static final String PHENODATA_ICON = "/eclipse/phenodata.gif";
	public static final String IMPORT_CANCEL_ICON = "/eclipse/import_cancel.gif";
	public static final String IMPORT_NEXT_ICON = "/eclipse/import_forward.gif";
	public static final String IMPORT_FINISH_ICON = "/eclipse/import_finish.png";
	public static final String IMPORT_BACK_ICON = "/eclipse/import_backward.gif";
	public static final String IMPORT_RESET_ICON = "/eclipse/import_reset.png";
	public static final String IMPORT_HEADER_ICON = "/eclipse/import_header.png";
	public static final String IMPORT_FOOTER_ICON = "/eclipse/import_footer.png";
	public static final String IMPORT_TITLE_ICON = "/eclipse/import_title.png";
	public static final String DEFAULT_VIEW_MENUICON = "/eclipse/defaultView.gif";
	public static final String EXPORT_MENUICON = "/eclipse/export.gif";
	public static final String DELETE_MENUICON = "/eclipse/delete.gif";
	public static final String HELP_MENUICON = "/eclipse/help.gif";
	public static final String LINK_PHENODATA_MENUICON = "/link_phenodata.png";
	public static final String LINK_MENUICON = "/eclipse/link.gif";
	public static final String UNLINK_MENUICON = "/eclipse/unlink.gif";
	public static final String QUICKLINK_ICON = "/eclipse/bulb.gif";

	// rest of the icons done by us (mostly by Mikko Koski and Petri Klemelä)
	public static String ZOOM_IN_ICON = "/zoom-in.png";
	public static String ZOOM_OUT_ICON = "/zoom-out.png";
	public static String MAGNIFIER_ICON = "/viewmag.png";
	public static String CLOSE_FILE_ICON = "/fileclose.png";
	public static final String APPLICATION_ICON = "/chipster_icon.png";
	public static final String LOGIN_BANNER = "/login_banner.png";

	public static String ZOOM_IN_CURSOR_IMAGE = "/zoomInCursor.png";
	public static String ZOOM_OUT_CURSOR_IMAGE = "/zoomOutCursor.png";
	public static String ROTATE_CURSOR_IMAGE = "/3dRotateCursor.png";
	public static String ROTATE_AND_ZOOM_CURSOR_IMAGE = "/rotatezoom.png";
	public static String ROTATE_IMAGE = "/3dRotate.png";
	public static String ARROW_ICON = "/arrow.png";
	public static String HAND_ICON = "/hand.png";
	public static final String DOUBLE_FORWARD_ICON = "/forward.png";
	public static final String DOUBLE_FORWARD_BW_ICON = "/forward_bw.png";
	public static final String SUITABLE_ICON = "/yes.png";
	public static final String INCOMPATIBLE_ICON = "/no.png";
	public static final String SUITABILITY_WARNING_ICON = "/maybe.png";
	public static String GENERATE_HISTORY_ICON = "/history.png";
	public static String TO_TOP_ICON = "/toTop.png";
	public static String TO_BOTTOM_ICON = "/toBottom.png";
	public static String REDRAW_ICON = "/redraw_bw.png";
	public static String MAXIMISE_ICON = "/maximise.png";
	public static String RESTORE_ICON = "/restore.png";
	public static String QUESTION_MARK_ICON = "/question_mark.png";
	public static String TO_WINDOW_ICON = "/to_window.png";
	public static String CLOSE_ICON = "/close.png";

	public static String XY_PLANE = "/XYPlane.png";
	public static String XZ_PLANE = "/XZPlane.png";
	public static String YZ_PLANE = "/YZPlane.png";

	// Crystal icons, LGPL, http://www.everaldo.com/crystal/
	public static String PHENODATA_MENUICON = "/listicons/crystal-LGPL/kate.png";
	public static String TEXT_MENUICON = "/listicons/crystal-LGPL/txt.png";
	public static String IMAGE_MENUICON = "/listicons/crystal-LGPL/no3D.png";
	public static String HTML_MENUICON = "/listicons/crystal-LGPL/agt_web.png";

	// Modified icons based on the Crystal icons, LGPL,
	// http://www.everaldo.com/crystal/
	public static String PDF_MENUICON = "/listicons/crystal-LGPL/modified/man.png";
	public static String EXT_BROWSER_MENUICON = "/listicons/crystal-LGPL/modified/ext-browser.png";
	public static String SPREADSHEET_MENUICON = "/listicons/crystal-LGPL/modified/spreadsheet.png";
	public static String ARRAY_MENUICON = "/listicons/crystal-LGPL/modified/array.png";

	// rest of the icons done by us (mostly by Mikko Koski and Petri Klemelä)
	public static String VENN_MENUICON = "/listicons/venndi-48.png";

	public static String SOM_MENUICON = "/listicons/som2.png";
	public static String SCATTER_MENUICON = "/listicons/scatter4.png";
	public static String PROFILE_MENUICON = "/listicons/profile4.png";
	public static String SCATTER3DPCA_MENUICON = "/listicons/pca.png";
	public static String HC_MENUICON = "/listicons/hc3.png";
	public static String SCATTER3D_MENUICON = "/listicons/scatter3d_2.png";
	public static String PROFILES_MENUICON = "/listicons/profiles3.png";
	public static String VOLCANO_MENUICON = "/listicons/volcano2.png";
	public static String HISTOGRAM_MENUICON = "/listicons/histogram2.png";
	public static String GB_MENUICON = "/listicons/gb3.png";
	public static String HEATMAP_MENUICON = "/listicons/heatmap.png";

	public static final String SPLASH_SCREEN = "/splash.png";
	public static final String SPLASH_SCREEN_KIELIPANKKI = "/mylly-splash.png";
	public static final Color SPLASH_BORDER_COLOR = new Color(150, 150, 150);
	public static final Font SPLASH_SCREEN_FONT = new Font("SansSerif",
			Font.PLAIN, 10);

	public static final Font VISUALISATION_TITLE_FONT = new Font("SansSerif",
			Font.PLAIN, 16);

	public static final Font MONOSPACED_FONT = new Font("Monospaced",
			Font.PLAIN, 11);

	public static final String ICON_TYPE_BINARY = "/types/binary.png";
	public static final String ICON_TYPE_FOLDER = "/types/folder.gif";
	public static final String ICON_TYPE_HTML = "/types/html.png";
	public static final String ICON_TYPE_IMAGE = "/types/image.png";
	public static final String ICON_TYPE_PHENODATA = "/types/phenodata.png";
	public static final String ICON_TYPE_RAWDATA = "/types/rawdata.png";
	public static final String ICON_TYPE_TABLE = "/types/table.png";
	public static final String ICON_TYPE_TEXT = "/types/text.png";

	public static final String ARROW_UP_ICON = "/arrow_up.png";
	public static final String ARROW_DOWN_ICON = "/arrow_down.png";
	public static final String STOP_ICON = "/stop.png";
	public static final String RUNNING_ICON = "/running.gif";

	public static final String OPEN_SESSION_LINK_ICON = "/tree.png";
	public static final String IMPORT_LINK_ICON = "/table.png";
	public static final String EXAMPLE_SESSION_ICON = "/try.png";
    

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

	public static ImageIcon getIcon(String iconPath) {
		if (iconPath == null) {
			return null;
		}
		return new ImageIcon(VisualConstants.class.getResource(iconPath));		
	}
}
