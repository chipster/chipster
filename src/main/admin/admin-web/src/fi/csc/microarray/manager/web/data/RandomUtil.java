package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Random;

public class RandomUtil {
	
	private static final String[] usernames = new String[] { "korhonen", "virtanen", "makinen", "nieminen", "makela", "hamalain", "laine", "heikkine", "koskinen", "jarvinen" };

	
    public static Date getRandomDate(Random rnd) {
    	Calendar cal = new GregorianCalendar();
    	
    	cal.set(rnd.nextInt(5) + 2007, rnd.nextInt(12) + 1, rnd.nextInt(28) + 1, rnd.nextInt(24), rnd.nextInt(60), rnd.nextInt(60));
    	
    	return cal.getTime();
    }
    
    public static String getRandomUserName(Random rnd) {
        return usernames[rnd.nextInt(usernames.length)];
        
    }

}
