package fi.csc.microarray.client.visualisation.methods.gbrowser.utils;
import java.io.File;
import java.io.RandomAccessFile;


public class FileReadSpeedTest {

	public void read(int chunkSize, long seek, long count){	

		try {
			RandomAccessFile raf = new RandomAccessFile(new File("2000.fsf"), "r");

			byte[] chunk = new byte[chunkSize];			

			long start = System.currentTimeMillis();
			int i = 0;

			while( i < count){

				raf.read(chunk);
				raf.seek(seek);						

				i++;
			}

			System.out.println(System.currentTimeMillis() - start + "\t" + chunkSize + "\t" + seek + "\t" + count);
			
			raf.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args){
		
		int chunk = 8*1024;
			for(long seek = 2; seek <= 1024*1024*1024; seek*=2){
				long count = 2*1000*1000*1000 / (seek + chunk);
				(new FileReadSpeedTest()).read(chunk, seek, count);
			}
		
	}
}
