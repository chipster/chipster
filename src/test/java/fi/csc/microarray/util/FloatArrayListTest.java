package fi.csc.microarray.util;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class FloatArrayListTest {

	@Test
	public void test() {
		int size = 100000;
		int iterations = 1000;
		int tuneDownDivider = 10; // otherwise ArrayList would be too slow
		float[] a = new float[size];		
		List<Float> b = new ArrayList<Float>(size);
		FloatArrayList c = new FloatArrayList(size); 
		for (int i = 0; i < size; i++) {
			b.add(0.0f);
			c.add(0.0f);			
		}
		
		long time1 = System.currentTimeMillis();
		for (int r = 0; r < iterations; r++) {
			for (int i = 0; i < a.length; i++) {
				a[i] = a[i] + 2.0f;
			}
		}
		float diff1 = (float)(System.currentTimeMillis()-time1);
		
		long time2 = System.currentTimeMillis();
		for (int r = 0; r < iterations/tuneDownDivider; r++) {
			for (int i = 0; i < b.size(); i++) {
				b.set(i, b.get(i) + 2.0f);
			}
		}
		float diff2 = (float)(System.currentTimeMillis()-time2)*tuneDownDivider;
		
		long time3 = System.currentTimeMillis();
		for (int r = 0; r < iterations; r++) {
			for (int i = 0; i < c.size(); i++) {
				c.setFloat(i, c.get(i) + 2.0f);
			}
		}
		float diff3 = (float)(System.currentTimeMillis()-time3);
		
		// verify that data is correct
		float result = 2.0f*(float)iterations;
		Assert.assertTrue(a.length == b.size());
		Assert.assertTrue(b.size() == c.size());
		for (int i = 0; i < a.length; i++) {
			Assert.assertTrue("got " + a[i], closeEnough(a[i], result));
			Assert.assertTrue("got " + b.get(i)*(float)tuneDownDivider, closeEnough(b.get(i)*(float)tuneDownDivider, result));
			Assert.assertTrue("got " + c.get(i), closeEnough(c.get(i), result));
		}
		
		Assert.assertTrue("FloatArrayList was " + diff3/diff1 + " times slower than primitive array, which is too much", diff3/diff1 < 15.0f);
		Assert.assertTrue("FloatArrayList was only " + diff2/diff3 + " times faster than ArrayList, which is too little", diff2/diff3 > 3.0f);
	}
	
	private boolean closeEnough(float f1, float f2) {
		return Math.abs(f1-f2) < 0.01f;
	}
}
