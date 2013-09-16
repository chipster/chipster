/*
 * Copyright (c) 2000-2006 JGoodies Karsten Lentzsch. All Rights Reserved.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *  o Redistributions of source code must retain the above copyright notice, 
 *    this list of conditions and the following disclaimer. 
 *     
 *  o Redistributions in binary form must reproduce the above copyright notice, 
 *    this list of conditions and the following disclaimer in the documentation 
 *    and/or other materials provided with the distribution. 
 *     
 *  o Neither the name of JGoodies Karsten Lentzsch nor the names of 
 *    its contributors may be used to endorse or promote products derived 
 *    from this software without specific prior written permission. 
 *     
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */

package com.jgoodies.uif_lite.panel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.LayoutManager;
import java.awt.Paint;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JToolBar;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.AbstractBorder;

import org.apache.log4j.Logger;

/** 
 * A <code>JPanel</code> subclass that has a drop shadow border and 
 * that provides a header with icon, title and tool bar.<p>
 * 
 * This class can be used to replace the <code>JInternalFrame</code>,
 * for use outside of a <code>JDesktopPane</code>. 
 * The <code>SimpleInternalFrame</code> is less flexible but often
 * more usable; it avoids overlapping windows and scales well 
 * up to IDE size.
 * Several customers have reported that they and their clients feel 
 * much better with both the appearance and the UI feel.<p>
 * 
 * The SimpleInternalFrame provides the following bound properties:
 * <i>frameIcon, title, toolBar, content, selected.</i><p>
 * 
 * By default the SimpleInternalFrame is in <i>selected</i> state.
 * If you don't do anything, multiple simple internal frames will
 * be displayed as selected.
 * 
 * @author Karsten Lentzsch
 * @version $Revision: 1.11 $
 * 
 * @see    javax.swing.JInternalFrame
 * @see    javax.swing.JDesktopPane
 */

public class SimpleInternalFrame extends JPanel implements MouseListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(SimpleInternalFrame.class);

	private final int SPLIT_MARGIN = 6;
	
    private JLabel        titleLabel;
    private JPanel		  gradientPanel;
    private JPanel        headerPanel;
    private boolean       selected;
    private boolean       paintGradient;
    private boolean 	  paintTitleBorder;
    
    
    // Instance Creation ****************************************************

    /**
     * Constructs a SimpleInternalFrame with the specified title.
     * The title is intended to be non-blank, or in other words
     * should contain non-space characters.
     * 
     * @param title       the initial title
     */
    public SimpleInternalFrame(String title) {
        this(null, title, null, null);
    }
    
    
    /**
     * Constructs a SimpleInternalFrame with the specified 
     * icon, and title.
     * 
     * @param icon        the initial icon
     * @param title       the initial title
     */
    public SimpleInternalFrame(Icon icon, String title) {
        this(icon, title, null, null);
    }

    
    /**
     * Constructs a SimpleInternalFrame with the specified 
     * title, tool bar, and content panel.
     * 
     * @param title       the initial title
     * @param bar         the initial tool bar
     * @param content     the initial content pane
     */
    public SimpleInternalFrame(String title, JToolBar bar, JComponent content) {
        this(null, title, bar, content);
    }
    

    /**
     * Constructs a SimpleInternalFrame with the specified 
     * icon, title, tool bar, and content panel.
     * 
     * @param icon        the initial icon
     * @param title       the initial title
     * @param bar         the initial tool bar
     * @param content     the initial content pane
     */
    public SimpleInternalFrame(
        Icon icon,
        String title,
        JToolBar bar,
        JComponent content) {
        super(new BorderLayout());
        this.selected = false;
        this.titleLabel = new JLabel(title, icon, SwingConstants.LEADING);        
        titleLabel.setFont(getTittleFont());        
        
        // Minimizing by double clicking
        this.addMouseListener(this);
        
        this.paintGradient = true;
        this.paintTitleBorder = false;
        
        JPanel top = buildHeader(titleLabel, bar);

        add(top, BorderLayout.NORTH);
        if (content != null) {
            setContent(content);
        }
        setBorder(new ShadowBorder());
        setSelected(true);
        updateHeader();
    }

    
    // Public API ***********************************************************

    /**
     * Returns the frame's icon.
     * 
     * @return the frame's icon
     */
    public Icon getFrameIcon() {
        return titleLabel.getIcon();
    }
    

    /**
     * Sets a new frame icon.
     * 
     * @param newIcon   the icon to be set
     */
    public void setFrameIcon(Icon newIcon) {
        Icon oldIcon = getFrameIcon();
        titleLabel.setIcon(newIcon);
        firePropertyChange("frameIcon", oldIcon, newIcon);
    }
    

    /**
     * Returns the frame's title text.     *      * @return String   the current title text
     */
    public String getTitle() {
        return titleLabel.getText();
    }
    
    
    /**
     * Sets a new title text.     *      * @param newText  the title text tp be set
     */
    public void setTitle(String newText) {
        String oldText = getTitle();
        titleLabel.setText(newText);        
        firePropertyChange("title", oldText, newText);
    }
    
    /**
     * Returns the current toolbar, null if none has been set before.
     * 
     * @return the current toolbar - if any
     */
    public JToolBar getToolBar() {
        return headerPanel.getComponentCount() > 1
            ? (JToolBar) headerPanel.getComponent(1)
            : null;
    }
    

    /**
     * Sets a new tool bar in the header.
     * 
     * @param newToolBar the tool bar to be set in the header
     */
    public void setToolBar(JToolBar newToolBar) {
        JToolBar oldToolBar = getToolBar();
        if (oldToolBar == newToolBar) {
            return;
        }
        if (oldToolBar != null) {
            headerPanel.remove(oldToolBar);
        }
        if (newToolBar != null) {
            newToolBar.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
            headerPanel.add(newToolBar, BorderLayout.SOUTH);
        }
        updateHeader();
        firePropertyChange("toolBar", oldToolBar, newToolBar);
    }

    
    /**
     * Returns the content - null, if none has been set.
     *      * @return the current content     */
      
    public Component getContent() {
        return hasContent() ? getComponent(1) : null;
    }
    
    
    /**
     * Sets a new panel content; replaces any existing content, if existing.
     *      * @param newContent   the panel's new content     */
    public void setContent(Component newContent) {
        Component oldContent = getContent();
        if (hasContent()) {
            remove(oldContent);
        }
        add(newContent, BorderLayout.CENTER);
        firePropertyChange("content", oldContent, newContent);
    }
    

    /**
     * Answers if the panel is currently selected (or in other words active)
     * or not. In the selected state, the header background will be
     * rendered differently.
     * 
     * @return boolean  a boolean, where true means the frame is selected 
     *                  (currently active) and false means it is not  
     */
    public boolean isSelected() {
        return selected;
    }
    
    /**
     * Returns true if the frame is maximized
     * @return
     */
    public boolean isMaximized(){
    	if(this.getParent() instanceof JSplitPane){
    		JSplitPane split = ((JSplitPane)this.getParent());
    		if(split.getTopComponent() == this){
    			if(split.getDividerLocation() == split.getHeight() - (gradientPanel.getHeight() + split.getDividerSize() + SPLIT_MARGIN)){
    				logger.debug("is a maximized top component");
    				return true;
    			} else {
    				logger.debug("is not a maximized component, but it is a top component");
    				return false;
    			}
    		} else {
    			if(split.getDividerLocation() == gradientPanel.getHeight()){
    				logger.debug("is a maximized bottom component");
    				return true;
    			} else {
    				logger.debug("is not a maximized component, but is is a bottom component");
    				return false;
    			}
    		}
    	} else {
    		throw new IllegalStateException("SimpleInternalFrame is not on a split pane");
    	}
    }
    
    
    /**
     * This panel draws its title bar differently if it is selected,
     * which may be used to indicate to the user that this panel
     * has the focus, or should get more attention than other
     * simple internal frames.
     *
     * @param newValue  a boolean, where true means the frame is selected 
     *                  (currently active) and false means it is not
     */
    public void setSelected(boolean newValue) {
        boolean oldValue = isSelected();
        selected = newValue;
        updateHeader();
        firePropertyChange("selected", oldValue, newValue);
    }
    

    // Building *************************************************************

    /**
     * Creates and answers the header panel, that consists of:
     * an icon, a title label, a tool bar, and a gradient background.
     * 
     * @param label   the label to paint the icon and text
     * @param bar     the panel's tool bar
     * @return the panel's built header area
     */
    private JPanel buildHeader(JLabel label, JToolBar bar) {
        if(paintGradient){
        	gradientPanel =
        		new GradientPanel(new BorderLayout(), getHeaderBackground());
        } else {
        	gradientPanel = new JPanel(new BorderLayout());
        }
        label.setOpaque(false);

        gradientPanel.add(label, BorderLayout.WEST);
        gradientPanel.setBorder(BorderFactory.createEmptyBorder(3, 4, 3, 1));
        headerPanel = new JPanel(new BorderLayout());
        headerPanel.add(gradientPanel, BorderLayout.CENTER);
        setToolBar(bar);
        if(paintTitleBorder){
        	headerPanel.setBorder(new RaisedHeaderBorder());
        }
        headerPanel.setOpaque(false);
        return headerPanel;
    }

    /**
     * Updates the header.
     */
    private void updateHeader() {
        gradientPanel.setBackground(getHeaderBackground());
        //gradientPanel.setOpaque(isSelected());
        titleLabel.setFont(getTittleFont());
        titleLabel.setForeground(getTextForeground(isSelected()));
        headerPanel.repaint();
    }
    

    /**
     * Updates the UI. In addition to the superclass behavior, we need
     * to update the header component.
     */
    public void updateUI() {
        super.updateUI();
        if (titleLabel != null) {
            updateHeader();
        }
    }


    // Helper Code **********************************************************

    /**
     * Checks and answers if the panel has a content component set.
     * 
     * @return true if the panel has a content, false if it's empty
     */
    private boolean hasContent() {
        return getComponentCount() > 1;
    }
        /**
     * Determines and answers the header's text foreground color.
     * Tries to lookup a special color from the L&amp;F.
     * In case it is absent, it uses the standard internal frame forground.
     * 
     * @param isSelected   true to lookup the active color, false for the inactive
     * @return the color of the foreground text
     */
    protected Color getTextForeground(boolean isSelected) {
        Color c =
            UIManager.getColor(
                isSelected
                    ? "SimpleInternalFrame.activeTitleForeground"
                    : "SimpleInternalFrame.inactiveTitleForeground");
        if (c != null) {
            return c;
        }
        return UIManager.getColor(
            isSelected 
                ? "InternalFrame.activeTitleForeground" 
                : "Label.foreground");

    }
    
    protected Font getTittleFont() {       
        return UIManager.getFont("Label.font").deriveFont(Font.BOLD);

    }

    /**
     * Determines and answers the header's background color.
     * Tries to lookup a special color from the L&amp;F.
     * In case it is absent, it uses the standard internal frame background.
     * 
     * @return the color of the header's background
     */
    protected static Color getHeaderBackground() {

        Color c =
            UIManager.getColor("SimpleInternalFrame.activeTitleBackground");
        return c != null
            ? c
            : UIManager.getColor("InternalFrame.activeTitleBackground");
    }


    // Helper Classes *******************************************************

    // A custom border for the raised header pseudo 3D effect.
    private static class RaisedHeaderBorder extends AbstractBorder {
    	
    	private static Insets INSETS = new Insets(1, 1, 1, 0);
    		
        public Insets getBorderInsets(Component c) { return INSETS; }

        public void paintBorder(Component c, Graphics g,
            int x, int y, int w, int h) {

        	g.translate(x, y);
        	g.setColor(UIManager.getColor("controlLtHighlight"));
        	g.fillRect(0, 0,   w, 1);
        	g.fillRect(0, 1,   1, h-1);
        	g.setColor(UIManager.getColor("controlShadow"));
        	g.fillRect(0, h-1, w, 1);
        	g.translate(-x, -y);        	
        }
    }

    // A custom border that has a shadow on the right and lower sides.
    private static class ShadowBorder extends AbstractBorder {

        private static final Insets INSETS = new Insets(1, 1, 3, 3);

        public Insets getBorderInsets(Component c) { return INSETS; }

        public void paintBorder(Component c, Graphics g,
            int x, int y, int w, int h) {
                
            Color shadow        = UIManager.getColor("controlShadow");
            if (shadow == null) {
                shadow = Color.GRAY;
            }
            Color lightShadow   = new Color(shadow.getRed(), 
                                            shadow.getGreen(), 
                                            shadow.getBlue(), 
                                            170);
            Color lighterShadow = new Color(shadow.getRed(),
                                            shadow.getGreen(),
                                            shadow.getBlue(),
                                            70);
            g.translate(x, y);
            
            g.setColor(shadow);
            g.fillRect(0, 0, w - 3, 1);
            g.fillRect(0, 0, 1, h - 3);
            g.fillRect(w - 3, 1, 1, h - 3);
            g.fillRect(1, h - 3, w - 3, 1);
            // Shadow line 1
            g.setColor(lightShadow);
            g.fillRect(w - 3, 0, 1, 1);
            g.fillRect(0, h - 3, 1, 1);
            g.fillRect(w - 2, 1, 1, h - 3);
            g.fillRect(1, h - 2, w - 3, 1);
            // Shadow line2
            g.setColor(lighterShadow);
            g.fillRect(w - 2, 0, 1, 1);
            g.fillRect(0, h - 2, 1, 1);
            g.fillRect(w-2, h-2, 1, 1);
            g.fillRect(w - 1, 1, 1, h - 2);
            g.fillRect(1, h - 1, w - 2, 1);
            g.translate(-x, -y);
        }
    }

    
    /**
     * A panel with a horizontal gradient background.
     */
    private static final class GradientPanel extends JPanel {
        
        private GradientPanel(LayoutManager lm, Color background) {
            super(lm);
            setBackground(background);
        }

        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (!isOpaque()) {
                return;
            }
            //Color control = UIManager.getColor("control");
            //Color end = new Color((getBackground().getRed() + control.getRed()) / 2, (getBackground().getGreen() + control.getGreen()) / 2, (getBackground().getBlue() + control.getBlue()) / 2);
            int width  = getWidth();
            int height = getHeight();

            Graphics2D g2 = (Graphics2D) g;
            Paint storedPaint = g2.getPaint();
          /*  g2.setPaint(
                new GradientPaint(0, 0, getBackground(), width, 0, end));
            	//new GradientPaint(0, height/2, getBackground(), 0, height, VisualConstants.TITLECOLOR_BLUE.darker()));            
            	//new GradientPaint(0, height/2, getBackground(), 0, height, getBackground().darker()));
            g2.fillRect(0, 0, width, height);*/
            

           multipleGradients(g2,width,height, getHeaderBackground());
                              
            g2.setPaint(storedPaint);
        }
    }
    
    
    public static void multipleGradients(Graphics2D g2, int width, int height, Color background){
    	 Color startColor;
         Color endColor;
         int startHeight;
         int endHeight;                
         
         startColor = brighten(brighten(brighten(background)));
         endColor = background;
         startHeight = (int)(height*0.0/4.0);
         endHeight = (int)(height*1.0/4.0);
                     
         g2.setPaint(new GradientPaint(0, startHeight, startColor, 0, endHeight, endColor));
         g2.fillRect(0, startHeight, width, endHeight);
         
         startColor = background;
         endColor = background;
         startHeight = (int)(height*1.0/4.0);
         endHeight = (int)(height*2.0/4.0);
                     
         g2.setPaint(new GradientPaint(0, startHeight, startColor, 0, endHeight, endColor));
         g2.fillRect(0, startHeight, width, endHeight);
         
         startColor = background;
         endColor = brighten(background);
         startHeight = (int)(height*2.0/4.0);
         endHeight = (int)(height*3.0/4.0);
                     
         g2.setPaint(new GradientPaint(0, startHeight, startColor, 0, endHeight, endColor));
         g2.fillRect(0, startHeight, width, endHeight);
         
         startColor = brighten(background);
         endColor = darken(darken(background));
         startHeight = (int)(height*3.0/4.0);
         endHeight = (int)(height*4.0/4.0);
                     
         g2.setPaint(new GradientPaint(0, startHeight, startColor, 0, endHeight, endColor));
         g2.fillRect(0, startHeight, width, endHeight);
    }
    
    private static final double FACTOR = 0.9;
    
    public static Color darken(Color c) {
    	return new Color(Math.max((int)(c.getRed()  *FACTOR), 0), 
    			Math.max((int)(c.getGreen()*FACTOR), 0),
    			Math.max((int)(c.getBlue() *FACTOR), 0));
    }
    
    public static Color brighten(Color c) {
        int r = c.getRed();
        int g = c.getGreen();
        int b = c.getBlue();

        int i = (int)(1.0/(1.0-FACTOR));
        if ( r == 0 && g == 0 && b == 0) {
           return new Color(i, i, i);
        }
        if ( r > 0 && r < i ) r = i;
        if ( g > 0 && g < i ) g = i;
        if ( b > 0 && b < i ) b = i;

        return new Color(Math.min((int)(r/FACTOR), 255),
                         Math.min((int)(g/FACTOR), 255),
                         Math.min((int)(b/FACTOR), 255));
    }


    /**
     * Maximizes the frame if it is on a split pane and user double clicks 
     * to the title panel
     */
	public void mouseClicked(MouseEvent e) {
		
		if(SwingUtilities.isLeftMouseButton(e) && e.getClickCount() == 2){
			if(this.getParent() instanceof JSplitPane){
				JSplitPane split = (JSplitPane)this.getParent();
				
				int maximizedDividerLocation = 0;
				
				if(split.getTopComponent() == this) {
					maximizedDividerLocation = split.getHeight() - (gradientPanel.getHeight() + split.getDividerSize() + SPLIT_MARGIN);
				} else {
					maximizedDividerLocation = gradientPanel.getHeight();
				}
				
				if (isMaximized()) {
					split.setDividerLocation(split.getLastDividerLocation());
				} else {
					split.setDividerLocation(maximizedDividerLocation);
				}
			}
		}
	}


	public void mouseEntered(MouseEvent e) {
		// Do nothing
	}


	public void mouseExited(MouseEvent e) {
		// Do nothing
	}


	public void mousePressed(MouseEvent e) {
		// Do nothing
	}


	public void mouseReleased(MouseEvent e) {
		// Do nothing
	}
}