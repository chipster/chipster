package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

public class Utils {
	
	public static String toHumanReadable(long i){
		return toHumanReadable(i, true);
	}

	public static String toHumanReadable(long i, boolean returnZero) {

		if(i == 0){
			if(returnZero){
				return "0";
			} else {
				return "";
			}
		} else if (i < 0){
			return "" + i;
		}

		int pow = (int)Math.log10(i);

		String sym = "";
		if(pow >= 3){
			sym = "k";
		} 
		if(pow >= 6){
			sym = "M";
		} 
		if(pow >= 9){
			sym = "G";
		} 
		if(pow >= 12){
			sym = "T";
		} 

		int div = (int)Math.pow(10,  (pow - pow % 3));

		return "" + i / div + sym + " " + toHumanReadable(i % div, false);
	}

}
