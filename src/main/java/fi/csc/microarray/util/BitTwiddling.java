package fi.csc.microarray.util;

import java.math.BigInteger;
import java.util.Set;

public class BitTwiddling {

	public static final int BYTES_IN_INT = Integer.SIZE/Byte.SIZE;
	public static final int MASK_TO_BYTE = 0xFF;
	
	public static int[] calculateBitFrequencies(Set<Long> numbers) {
		int[] bitCounts = new int[Long.SIZE];		
		for (long value : numbers) {
			for (int b = 0; b < bitCounts.length; b++) {
				if (BigInteger.valueOf(value).testBit(b)) {
					bitCounts[b]++;
				}
			}
		}
		return bitCounts;
	}	
	
	
	public static byte[] intToBytes(int value) {
		
		byte[] bytes = new byte[4];
	   
	    // unrolled loop of 4 iterations
	    bytes[0] = (byte)(value & MASK_TO_BYTE);
	    value >>= Byte.SIZE;

	    bytes[1] = (byte)(value & MASK_TO_BYTE);
	    value >>= Byte.SIZE;

	    bytes[2] = (byte)(value & MASK_TO_BYTE);
	    value >>= Byte.SIZE;

	    bytes[3] = (byte)(value & MASK_TO_BYTE);
	    
	    return bytes;
	}

	   
	public static int bytesToInt(byte[] bytes)
	{
	    assert bytes.length >= BYTES_IN_INT : "bytes should have size " + BYTES_IN_INT + ", not " + bytes.length;

	    // unrolled loop of 4 iterations
	    int result = (bytes[0] & MASK_TO_BYTE);
	    result |= ((bytes[1] & MASK_TO_BYTE) << 8);
	    result |= ((bytes[2] & MASK_TO_BYTE) << 16);
	    result |= ((bytes[3] & MASK_TO_BYTE) << 24);

	    return result;
	}

}
