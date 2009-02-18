package fi.csc.microarray.util;

import java.math.BigInteger;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.security.SecureSessionPool;
import fi.csc.microarray.security.SecureSessionPool.Session;

public class SecureSessionPoolTest  {
	
	@Test(groups = {"unit"} )
	public void testBasicUsage() {
		SecureSessionPool ssp = new SecureSessionPool();
		Session session = ssp.createSession();
		session.putParameter("one", 1);
		session.putParameter("two", 2);
		String id = session.getID();
		// we could pass id as a String over network etc. 
		Session reloadedSession = ssp.getSession(id);
		Assert.assertTrue(((Integer)reloadedSession.getParameter("one")) == 1);
		Assert.assertTrue(((Integer)reloadedSession.getParameter("two")) == 2);
		ssp.removeSession(reloadedSession);
		Assert.assertTrue(ssp.getSession(reloadedSession.getID()) == null);
	}

	@Test(groups = {"unit"} )
	public void testKeyRandomness() {
		SecureSessionPool ssp = new SecureSessionPool();

		// create random keys with SessionPool
		final int valCount = 1000;
		Set<Long> numbers = new HashSet<Long>();
		for (int i = 0; i < valCount; i++) {
			UUID id = UUID.fromString(ssp.createSession().getID());
			long xorredValue = id.getLeastSignificantBits() ^ id.getMostSignificantBits();
			Assert.assertFalse(numbers.contains(xorredValue), "duplicate key generated");
			numbers.add(xorredValue);
		}
		
		
		// calculate how many times different bits are up
		int[] bitCounts = BitTwiddling.calculateBitFrequencies(numbers);

		// check that all bits vary about as much
		final long average = valCount/2;
		final long allowedError = valCount/10; // 10% error is ok
		final long allowedErrorInAvg = valCount/100; // 1% error is ok
		long sum = 0;
		for (int i = 0; i < bitCounts.length; i++) {
			int count = bitCounts[i];
			sum += count;
			Assert.assertTrue(BigInteger.valueOf(count - average).abs().intValue() < allowedError, "bit count too far from theoretical average: " + count);
		}
		Assert.assertTrue(BigInteger.valueOf(sum/Long.SIZE-average).abs().intValue() < allowedErrorInAvg, "average too far from theoretical average");		
	}
	
}
