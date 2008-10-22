package fi.csc.microarray.util;

import java.awt.Color;

public class Colors {
	
	public static int HUE_COMPONENT = 0;
	public static int SATURATION_COMPONENT = 1;
	public static int BRIGHTNESS_COMPONENT = 2;
	
	public static Color adjust(Color original, float multiplier, int component) {		
		int r = original.getRed();
		int g = original.getGreen();
		int b = original.getBlue();
		float[] hsb = Color.RGBtoHSB(r, g, b, null);
		hsb[component] = multiplier*hsb[component]; // adjust 
		return new Color(Color.HSBtoRGB(hsb[0], hsb[1], hsb[2]));
	}
}
