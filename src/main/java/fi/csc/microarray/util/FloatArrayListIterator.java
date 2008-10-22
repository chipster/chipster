package fi.csc.microarray.util;

import java.util.Iterator;

public interface FloatArrayListIterator extends Iterator<Float> {

		public boolean hasNext();
		public Float next();
		public float nextFloat();
		public void remove();
}