package fi.csc.microarray.client.visualisation.methods.threed;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.LinkedList;
import java.util.concurrent.PriorityBlockingQueue;

import javax.swing.JComponent;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.event.MouseInputListener;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.constants.VisualConstants;

/**
 * The class that takes care of the user manipulation of the scatterplot, basically
 * handles selections, rotations etc. The base of the implemention is Threed-visualisation
 * from Viski library (http://www.cs.helsinki.fi/group/viski/ ).
 * 
 * @author Petri Klemelä
 */
public class CoordinateArea extends JComponent 
implements ActionListener, MouseInputListener, MouseWheelListener, PropertyChangeListener  {
	
	private static final double KINETIC_SPEED_FACTOR = 0.002;
	
	private static long MOVE_TIME_LIMIT = 10; //ms

	private JMenuItem hideSelected;
    private JMenuItem showAll;
    private JMenuItem invertSelection;   
        
	public enum PaintMode {
		SPHERE("Sphere"), RECT("Rectangle");
		
		private String name;

		PaintMode(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	};
    
    public PaintMode paintMode = PaintMode.SPHERE;
    
    //Just a little toy
    private boolean kineticMoveMode = false;
    
    public PaintMode getPaintMode() {
    	return paintMode;
    }
    
    public void setPaintMode(PaintMode mode) {
    	paintMode = mode;
    	this.repaint();
    }    

	private LinkedList<DataPoint> selectedPoints;
	private Projection projection;
	private Worker worker;
	protected AutomatedMovement movement;
	
	private ClientApplication application = Session.getSession().getApplication();

	Scatterplot3D controller;

	//Any number to make rotation speed slow enough
	final double ANGLE_INCREMENT = Math.PI/(360*4);
	
	private AutomatedMovement.RotationTask kineticMovement;
	private long lastDragEventTime = -1;

	/**
	 * 
	 * @param controller 
	 */
	public CoordinateArea(Scatterplot3D controller) {
		this.controller = controller;
		this.createPopupMenu();
		
		addMouseListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(this);
		setOpaque(true);

		setFocusable(true);
		
		this.setBackground(Color.black);

		selectedPoints = new LinkedList<DataPoint>();
		projection = new Projection(controller.getDataModel());
		worker = new Worker(this, projection);
		worker.start();
		worker.workRequest();
				
		movement = new AutomatedMovement(projection, worker);
		movement.start();
		kineticMovement = movement.startKineticMove(25, 0.95);		
		
		this.updateSelectedFromApplication();
		application.addClientEventListener(this);
		
	}

	public Projection getProjection(){
		return projection;
	}
	
	/**
	 * 
	 * @param g 
	 */
	protected void paintComponent(Graphics g) {
//		long startTime = System.currentTimeMillis();
		
		if(this.getHeight() <= this.getWidth()){
			projection.setViewWindowWidth(projection.getViewWindowHeight() * 
				this.getWidth() / (double)this.getHeight());
		} else {
			projection.setViewWindowHeight(projection.getViewWindowWidth() * 
					this.getHeight() / (double)this.getWidth());
		}
		
		Graphics2D g2d = (Graphics2D)g;
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);							

		PriorityBlockingQueue<Drawable> points = projection.getResultPoints();				
		
		if (points != null && !points.isEmpty()) {
			
			g2d.setColor(this.getBackground());
			g2d.fillRect(0, 0, getWidth(), getHeight());
			
			Drawable p = null;
			while ((p = points.poll()) != null) {
				p.draw(g2d, getWidth(), getHeight(), getPaintMode());
			}

			if(this.mouseDragged){
				g2d.setColor(Color.WHITE);
				g2d.setStroke(VisualConstants.dashLine);
				int x = mousePressX < mouseX ? mousePressX : mouseX;
				int y = mousePressY < mouseY ? mousePressY : mouseY;
				int w = Math.abs(mouseX - mousePressX);
				int h = Math.abs(mouseY - mousePressY);
				
				g2d.setColor(this.getForeground());
				g2d.drawRect(x,y,w,h);
			}
		} else {
			worker.workRequest();
		}
		
//		long endTime = System.currentTimeMillis();			
	}

	//Methods required by the MouseInputListener interface.
	/**
	 * 
	 * @param e 
	 */
	public void mouseClicked(MouseEvent e) {

		//double deg = 0;
		switch (e.getButton()) {
		case MouseEvent.BUTTON1:
			if (e.isControlDown()) {
				addToSelections(e);
			}
			else {
				selectOne(e);
			}
			break;
		case MouseEvent.BUTTON2:
			break;
		case MouseEvent.BUTTON3:
			break;
		default:
		}

		worker.workRequest();
	}

	public void rotateWithDrag(int mouseX, int mouseY, MouseEvent e){		

		double xfactor = Math.abs(mouseX- mousePressX);
		double yfactor = Math.abs(mouseY- mousePressY);
		
		controller.stopAutoRotation();
		
		if(e.isShiftDown()){
			double deg = 
				projection.getZAxisRotation();
			if (mousePressX - mouseX < 0) {
				projection.setZAxisRotation(deg +ANGLE_INCREMENT*xfactor);
			} else if (mousePressX - mouseX > 0) {
				projection.setZAxisRotation(deg - ANGLE_INCREMENT*xfactor);
			}
			this.zoom((mouseY - mousePressY) / 100.0 );
		} else {
			if(!kineticMoveMode){
				kineticMovement.stop();
				
				double deg = 
					projection.getYAxisRotation();
				if (mousePressX - mouseX < 0) {
					projection.setYAxisRotation(deg + ANGLE_INCREMENT*xfactor);
				} else if (mousePressX - mouseX > 0) {
					projection.setYAxisRotation(deg -ANGLE_INCREMENT*xfactor);
				}

				deg = 
					projection.getXAxisRotation();
				if (mousePressY - mouseY < 0) {
					projection.setXAxisRotation(deg + ANGLE_INCREMENT*yfactor);
				} else if (mousePressY - mouseY > 0) {
					projection.setXAxisRotation(deg -ANGLE_INCREMENT*yfactor);	
				}
				kineticMovement.setAngleIncs(KINETIC_SPEED_FACTOR*(mouseY - mousePressY), KINETIC_SPEED_FACTOR*( mouseX-mousePressX), 0);
			}
		}

		mousePressX = mouseX;
		mousePressY = mouseY;
		
		lastDragEventTime = System.currentTimeMillis();

		if(!kineticMoveMode){
			worker.workRequest();
		}
	}
	

	public void moveWithDrag(int mouseX, int mouseY, MouseEvent e){        	

		double[] orig = projection.getPointOfView();
		
		//TODO
		//Folloving equation is a flight of fantasy and doesn't work even right
		double divider = //Math.abs(orig[2]) *
		projection.getDistanceOfProjectionPlaneFromOrigin() *
		projection.getDistanceOfProjectionPlaneFromOrigin();
		//projection.getViewWindowWidth();

		orig[0] -= (mouseX - mousePressX) / divider;
		orig[1] -= (mouseY - mousePressY) / divider;

		projection.setPointOfView(orig);

		mousePressX = mouseX;
		mousePressY = mouseY;

		worker.workRequest();
	}

	public void mouseMoved(MouseEvent e) { }
	public void mouseExited(MouseEvent e) { }

	/**
	 * 
	 * @param e 
	 */
	public void mouseEntered(MouseEvent e) { 
		this.requestFocus();
	}

	private int mousePressX = 0;
	private int mousePressY = 0;

	//For selection rectangle drawing
	private int mouseX;
	private int mouseY;

	/**
	 * 
	 * @param e 
	 */
	public void mousePressed(MouseEvent e) {
		this.requestFocus();
		
		kineticMovement.stop();
		
		mousePressX = e.getX();
		mousePressY = e.getY();
		
		//kineticMoveMode = false;
	}

	/**
	 * 
	 * @param e 
	 */
	public void mouseDragged(MouseEvent e) {
		if(controller.getTool() == Scatterplot3D.Tool.ROTATE){
			this.rotateWithDrag(e.getX(), e.getY(), e);
		} else if(controller.getTool() == Scatterplot3D.Tool.MOVE){
			this.moveWithDrag(e.getX(), e.getY(), e);
		} else if(controller.getTool() == Scatterplot3D.Tool.SELECT){
			mouseDragged = true;
			mouseX = e.getX();
			mouseY = e.getY();
			worker.workRequest();
		}
	}        

	boolean mouseDragged = false;

	/**
	 * 
	 * @param e 
	 */
	public void mouseReleased(MouseEvent e) {
		if(System.currentTimeMillis() - lastDragEventTime < MOVE_TIME_LIMIT){			
			kineticMovement = movement.restartKineticMove();			
		} else {
			kineticMovement.stop();
		}
		
		if (e.getButton() == MouseEvent.BUTTON1 && mouseDragged == true) {
			selectGroup(e, mousePressX, mousePressY);
			mouseDragged = false;
		
			
			//repaint();
			worker.workRequest();
		}
	}
	/**
	 * 
	 * @param e 
	 */
	public void mouseWheelMoved(MouseWheelEvent e) {
		this.zoom(e.getWheelRotation() / 10.0);
	}
	
	private void zoom(double value){
		double[] pov = projection.getPointOfView();
		pov[2] -= value;
		projection.setPointOfView(pov);	
		worker.workRequest();
	}


	private void clearSelections() {
		if (selectedPoints == null)
			return;

		for (DataPoint dp : selectedPoints) {
			if (dp != null)
				dp.selected = false;
		}
		selectedPoints.clear();
		controller.getAnnotateList().setSelectedListContentAsDataPoints(selectedPoints, this, false, controller.getFrame().getDatas().get(0));
	}

	private void selectOne(MouseEvent e) {
		clearSelections();
		selectedPoints =
			DataPoint.getNearest(e.getX(), e.getY(), controller.getDataModel().getDataArray(), 8);
		for (DataPoint dp : selectedPoints) {
			if (dp != null)
				dp.selected = true;
		}
		
		controller.getAnnotateList().setSelectedListContentAsDataPoints(selectedPoints, this, true, controller.getFrame().getDatas().get(0));
	}

	private void addToSelections(MouseEvent e) {
		LinkedList<DataPoint> newSelection =
			DataPoint.getNearest(e.getX(), e.getY(), controller.getDataModel().getDataArray(), 4);
		for (DataPoint dp : newSelection) {
			if (dp != null) {
				if (dp.selected == true) {
					dp.selected = false;
					selectedPoints.remove(dp);
				}
				else {
					dp.selected = true;
					selectedPoints.add(dp);
				}
			}
		}
		controller.getAnnotateList().setSelectedListContentAsDataPoints(selectedPoints, this, true, controller.getFrame().getDatas().get(0));
	}

	private void selectGroup(MouseEvent e, int x1, int y1) {
		clearSelections();
		selectedPoints =
			DataPoint.getGroup(x1, y1, 
					e.getX(), e.getY(), controller.getDataModel().getDataArray());
		for (DataPoint dp : selectedPoints) {
			if (dp != null)
				dp.selected = true;
		}
		controller.getAnnotateList().setSelectedListContentAsDataPoints(selectedPoints, this, true, controller.getFrame().getDatas().get(0));
	}
	
	class PopupListener extends MouseAdapter {
	    JPopupMenu popup;
	    
	    PopupListener(JPopupMenu popupMenu) {
	        popup = popupMenu;
	    }
	    
	    public void mousePressed(MouseEvent e) {
	        maybeShowPopup(e);
	    }
	    
	    public void mouseReleased(MouseEvent e) {
	        maybeShowPopup(e);
	    }
	    
	    private void maybeShowPopup(MouseEvent e) {
	        if (e.isPopupTrigger()) {
	            int x = e.getX();
	            int y = e.getY();
	            
	            if (selectedPoints.isEmpty()) {
	                hideSelected.setEnabled(false);
	                invertSelection.setEnabled(false);
	                showAll.setEnabled(false);
	            }
	            else {
	                hideSelected.setEnabled(true);
	                invertSelection.setEnabled(true);
	                showAll.setEnabled(true);
	            }
	            popup.show(e.getComponent(), x, y);
	        }
	    }
	}

	/**
	 * 
	 * @param e 
	 */
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == hideSelected) {
	        for (Drawable d : selectedPoints) {
	            d.hidden = true;
	        }
	    }
	    else if (e.getSource() == showAll) {
	        for (Drawable d : controller.getDataModel().getDataArray()) {
	            d.hidden = false;
	        }
	    }
	    else if (e.getSource() == invertSelection) {
	        selectedPoints = new LinkedList<DataPoint>();
	        for (Drawable d : controller.getDataModel().getDataArray()) {
	            d.selected = !d.selected;
	            if (d instanceof DataPoint && d.selected)
	                selectedPoints.add((DataPoint)d);
	        }
	    }
	    //repaint();
	    worker.workRequest();
	}

	public void createPopupMenu() {
	    
	    //Create the popup menu.
	    JPopupMenu popup = new JPopupMenu();

	    hideSelected = new JMenuItem("Hide selected");
	    hideSelected.addActionListener(this);
	    popup.add(hideSelected);
	    showAll = new JMenuItem("Show all points");
	    showAll.addActionListener(this);
	    popup.add(showAll);
	    invertSelection = new JMenuItem("InvertSelection");
	    invertSelection.addActionListener(this);
	    popup.add(invertSelection);
	    
	    
	    //Add listener to the text area so the popup menu can come up.
	    MouseListener popupListener = new PopupListener(popup);
	    this.addMouseListener(popupListener);
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if(evt instanceof SelectionEvent && !evt.getSource().equals(this) && 
				((SelectionEvent)evt).getData() == controller.getFrame().getDatas().get(0)){			
			updateSelectedFromApplication();
			this.repaint();
		}
	}

	private void updateSelectedFromApplication() {
		Drawable[] drawables = controller.getDataModel().getDataArray();
		this.clearSelections();
		for (int index : application.getSelectionManager().
				getSelectionManager(controller.getFrame().getDatas().get(0)).getSelectionAsRows()){
			
			for(Drawable drawable : drawables){
				if(drawable instanceof DataPoint){
					DataPoint dp = (DataPoint)drawable;
					if(dp.getIndex() == index){
						dp.selected = true;
						selectedPoints.add(dp);
						break;
					}
				}
			}
		}
		
		controller.getAnnotateList().setSelectedListContentAsDataPoints(selectedPoints, this, false, controller.getFrame().getDatas().get(0));
	}
}
