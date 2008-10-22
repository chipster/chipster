package fi.csc.microarray.proto;

import java.net.URL;

import javax.imageio.ImageIO;
import javax.imageio.stream.ImageInputStream;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
 
import com.sun.imageio.plugins.gif.GIFImageMetadata;
import com.sun.imageio.plugins.gif.GIFImageReader;
 
public class AnimatedGifTest {
    public static void main(String[] args) {
        try {
            JFrame frame = new JFrame();
            frame.setDefaultCloseOperation(JDialog.EXIT_ON_CLOSE);
            URL url = new URL("file:examples/animated.gif");
            frame.add(new JLabel(new ImageIcon(url)));
            frame.pack();
            frame.setVisible(true);
            // Read the delay time for the first frame in the animated gif
            GIFImageReader reader = (GIFImageReader) ImageIO.getImageReadersByFormatName("GIF").next();
            ImageInputStream iis = ImageIO.createImageInputStream(url.openStream());
            reader.setInput(iis);
            GIFImageMetadata md = (GIFImageMetadata) reader.getImageMetadata(0);            
            System.out.println("Time Delay = " + md.delayTime);
            reader.dispose();
            iis.close();
        } catch (Exception e) {e.printStackTrace();}
    }    
}
