package fi.csc.microarray.util;

import java.util.UUID;

import org.testng.Assert;
import org.testng.annotations.Test;

public class BitTwiddlingTest {
	
	@Test(groups = {"unit"} )
	public void testByteConversions() {
		// int -> byte
		for (int i = 0; i < 1000; i++) {
			int original = (int)UUID.randomUUID().getLeastSignificantBits();
			int converted = BitTwiddling.bytesToInt(BitTwiddling.intToBytes(original));
			Assert.assertTrue(original == converted, "comparing 0x" + Integer.toHexString(original) + " and 0x" + Integer.toHexString(converted));
		}
		
		// byte -> int
		for (int i = 0; i < 1000; i++) {
			byte[] original = new byte[BitTwiddling.BYTES_IN_INT];
			for (int c = 0; c < original.length; c++) {
				original[c] = (byte)UUID.randomUUID().getLeastSignificantBits();
			}
			byte[] converted = BitTwiddling.intToBytes(BitTwiddling.bytesToInt(original));
			for (int c = 0; c < original.length; c++) {
				Assert.assertTrue(original[c] == converted[c], "comparing " + original[c] + " and " + converted[c]);	
			}
		}
	}
	
	public static void main(String[] args) {
		//int i = 0x11223344;
		int i = 0xe49b24f3;
		byte[] b = BitTwiddling.intToBytes(i);
		int i2 = BitTwiddling.bytesToInt(b);
		System.out.println(Integer.toHexString(i) + " and " + Integer.toHexString(i2));
		
	}
}
