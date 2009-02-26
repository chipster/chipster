package fi.csc.microarray.util;

/////////////////////////////////////////////////////////
//  Bare Bones Browser Launch                          //
//  Version 1.5 (December 10, 2005)                    //
//  By Dem Pilafian                                    //
//  (modified by Aleksi Kallio)                        //
//  Supports: Mac OS X, GNU/Linux, Unix, Windows XP    //
//  Example Usage:                                     //
//     String url = "http://www.centerkey.com/";       //
//     BareBonesBrowserLaunch.openURL(url);            //
//  Public Domain Software -- Free to Use as You Like  //
/////////////////////////////////////////////////////////

import java.lang.reflect.Method;

public class BrowserLauncher {

	private static final String[] UNIX_BROWSERS = {"firefox", "seamonkey", "mozilla", "konqueror", "opera", "epiphany", "netscape" };

	public static void openURL(String url) throws Exception {
		String osName = System.getProperty("os.name");

		if (osName.startsWith("Mac OS")) {
			Class<?> fileMgr = Class.forName("com.apple.eio.FileManager");
			Method openURL = fileMgr.getDeclaredMethod("openURL", new Class[] {String.class});
			openURL.invoke(null, new Object[] {url});

		} else if (osName.startsWith("Windows")) {
			Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);

		} else { 
			
			// assume *nix			
			String availableBrowser = null;
			for (String browser : UNIX_BROWSERS) {
				if (Runtime.getRuntime().exec(new String[] {"which", browser}).waitFor() == 0) { 
					availableBrowser = browser;
					break;
				}
			}
			
			if (availableBrowser == null) {
				throw new Exception("Could not find web browser");
				
			} else {
				Runtime.getRuntime().exec(new String[] {availableBrowser, url});
			}
		}

	}

}
