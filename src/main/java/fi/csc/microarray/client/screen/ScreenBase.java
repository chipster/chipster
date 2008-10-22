/*
 * Created on Feb 14, 2005
 *
 */
package fi.csc.microarray.client.screen;


/**
 * @author akallio
 *
 */
public abstract class ScreenBase implements Screen {

	protected Object childScreenParameter;
	
	public void setChildScreenParameter(Object childScreenParameter) {
		this.childScreenParameter = childScreenParameter;
		update();
	}
	
	public void update() {
		// empty implementation
	}

}
