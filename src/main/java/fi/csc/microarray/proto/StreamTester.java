/*
 * Created on May 2, 2005
 *
 */
package fi.csc.microarray.proto;


/**
 * @author akallio
 *
 */
public class StreamTester {
//	/**
//	 * Logger for this class
//	 */
//	private static final Logger logger = Logger.getLogger(StreamTester.class);
//
//	public static int PACKET_SIZE = 1024;
//	public static int PACKET_COUNT = 16*1024; // 16 megs
//	public static char PACKET_CONTENT = 'a';
//	
//	public static void main(String[] args) {
//		if (args.length == 0) {
//			args = new String[] {""};
//		}
//		
//		System.out.println("testing streaming...");
//		
//		if (!"send".equals(args[0])) {
//			System.out.println("receiving...");
//			new Thread(new Runnable() {
//				public void run() {
//					receive();
//				}
//			}).start();
//			
//			// wait for the receiver to start, so it doesn't miss any messages
//			try {
//				Thread.sleep(2000);
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//			
//		}
//		
//		if (!"receive".equals(args[0])) {
//			System.out.println("sending...");
//			new Thread(new Runnable() {
//				public void run() {
//					send();
//				}
//			}).start();
//		}
//	}
//	
//	public static void receive() {
//		try {
//			MessagingEndpoint endpoint = new MessagingEndpoint(new NodeBase() {
//				public String getName() {
//					return "_test";
//				}    		
//			});
//			MessagingTopic topic  = endpoint.createTopic(Topics.Name.TEST2_TOPIC);
//			
//			InputStream in = topic.openJMSInputStream();
//			
//			byte[] input = new byte[PACKET_SIZE];
//			in.read(input);
//			long start = System.currentTimeMillis();
//			for (long i = 0; i < PACKET_COUNT-1; i++) {
//				logger.debug("pre-read");
//				in.read(input);
//				logger.debug("post-read");
//				if (input[0] != PACKET_CONTENT) {
//					throw new Exception("error in data, got " + input[0]);
//				}
//			}
//			
//			long time = (System.currentTimeMillis() - start);
//			float speed = ((float)PACKET_COUNT*PACKET_SIZE*8) / ((float)time / 1000) / (1024F*1024F);
//			
//			System.out.println("receiving took " + time + " millis, " + speed + " MB/s");
//			
//			in.close();
//			
//			endpoint.close();
//			System.out.flush();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//	}
//
//	public static void send() {
//		try {
//			MessagingEndpoint endpoint = new MessagingEndpoint(new NodeBase() {
//				public String getName() {
//					return "_test";
//				}    		
//			});
//			MessagingTopic topic  = endpoint.createTopic(Topics.Name.TEST2_TOPIC);
//			
//			OutputStream out = topic.openJMSOutputStream();
//			
//			byte[] test = new byte[PACKET_SIZE];
//			Arrays.fill(test, (byte)PACKET_CONTENT);
//			long start = System.currentTimeMillis();
//			for (long i = 0; i < PACKET_COUNT; i++) {
//				logger.debug("pre-write");
//				out.write(test);
//				logger.debug("post-write");
//			}
//			System.out.println("sending took " + (System.currentTimeMillis() - start) + " millis");
//			
//			out.close();
//			endpoint.close();
//			
//			System.out.flush();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//	}
//	
}
