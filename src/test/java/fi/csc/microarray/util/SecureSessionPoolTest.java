package fi.csc.microarray.util;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.security.SecureSessionPool;
import fi.csc.microarray.security.SecureSessionPool.Session;

public class SecureSessionPoolTest  {
	
	@Before
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"auth"}));
	}
	
	@Test
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

	@Test
	public void testKeyRandomness() {
		SecureSessionPool ssp = new SecureSessionPool();

		// create random keys with SessionPool
		final int valCount = 1000;
		Set<Long> numbers = new HashSet<Long>();
		for (int i = 0; i < valCount; i++) {
			UUID id = UUID.fromString(ssp.createSession().getID());
			long xorredValue = id.getLeastSignificantBits() ^ id.getMostSignificantBits();
			Assert.assertFalse("duplicate key generated", numbers.contains(xorredValue));
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
			Assert.assertTrue("bit count too far from theoretical average: " + count, BigInteger.valueOf(count - average).abs().intValue() < allowedError);
		}
		Assert.assertTrue("average too far from theoretical average", BigInteger.valueOf(sum/Long.SIZE-average).abs().intValue() < allowedErrorInAvg);		
	}
	
}
