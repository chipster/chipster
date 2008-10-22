/*
 * Created on Feb 23, 2005
 *
 */
package fi.csc.microarray.util;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;


/**
 * @author akallio
 */
public class FloatArrayList implements List<Float> {
	
	private static final int DEFAULT_SIZE = 1024;
	private static final int MINIMUM_SIZE = 32;

	private float[] data;
	private int lastElement;

	public class InternalFloatArrayListIterator implements FloatArrayListIterator {

		int index = 0;

		public boolean hasNext() {
			return index < lastElement; 
		}

		public Float next() {
			return data[index++];
		}

		public float nextFloat() {
			return data[index++];
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
		
	}

	public FloatArrayList() {
		this(DEFAULT_SIZE);
	}
	
	public FloatArrayList(int size) {
		if (size < MINIMUM_SIZE) {
			size = MINIMUM_SIZE; // size must never be 0
		}
		data = new float[size];
		lastElement = 0;
	}

	public FloatArrayList(float[] data) {
		this.data = data;
		lastElement = data.length;
	}
	
	public FloatArrayList(double[] data) {
		this(toFloatArray(data));
	}

	public FloatArrayList(Float[] data) {
		this(toFloatArray(data));
	}

	public FloatArrayList(Double[] data) {
		this(toFloatArray(data));
	}
	
	public FloatArrayList(List<Float> values) {
		this(values.toArray(new Float[0]));
	}

	/**
	 * We need this helper method because in constructors call to another 
	 * constructor must be the first stament. 
	 */
	private static float[] toFloatArray(double[] data) {
		float[] f = new float[data.length];
		for (int i = 0; i < data.length; i++) {
			f[i] = (float)data[i];
		}
		return f;
	}

	/**
	 * We need this helper method because in constructors call to another 
	 * constructor must be the first stament. 
	 */
	private static float[] toFloatArray(Float[] data) {
		float[] f = new float[data.length];
		for (int i = 0; i < data.length; i++) {
			f[i] = data[i].floatValue(); // lets be explicit, we could use autoboxing
		}
		return f;
	}

	/**
	 * We need this helper method because in constructors call to another 
	 * constructor must be the first stament. 
	 */
	private static float[] toFloatArray(Double[] data) {
		float[] f = new float[data.length];
		for (int i = 0; i < data.length; i++) {
			f[i] = data[i].floatValue();
		}
		return f;
	}

	private void internalAddElement(float element) {
		if (lastElement < data.length) {
			data[lastElement] = element;
			lastElement++;			
		} else {
			// we have to grow
			float[] old = data;
			data = new float[old.length*2];
			System.arraycopy(old, 0, data, 0, old.length);
			
			//	this should work without growing now
			internalAddElement(element); 
		}
	}
	
	public int size() {
		return this.lastElement;
	}

	public boolean isEmpty() {
		throw new UnsupportedOperationException();
	}

	public boolean contains(Object o) {
		throw new UnsupportedOperationException();
	}

	public FloatArrayListIterator floatIterator() {
		return new InternalFloatArrayListIterator();
	}
	
	public Iterator<Float> iterator() {
		return new InternalFloatArrayListIterator();
	}

	public Object[] toArray() {
		throw new UnsupportedOperationException();
	}

	public <T> T[] toArray(T[] a) {
		throw new UnsupportedOperationException();
	}

	public boolean add(Float element) {
		internalAddElement(element.floatValue());
		return true;
	}

	public boolean add(float element) {
		internalAddElement(element);
		return true;
	}

	public boolean remove(Object o) {
		throw new UnsupportedOperationException();
	}

	public boolean containsAll(Collection< ? > c) {
		throw new UnsupportedOperationException();
	}

	public boolean addAll(Collection< ? extends Float> c) {
		throw new UnsupportedOperationException();
	}

	public boolean addAll(int index, Collection< ? extends Float> c) {
		throw new UnsupportedOperationException();
	}

	public boolean removeAll(Collection< ? > c) {
		throw new UnsupportedOperationException();
	}

	public boolean retainAll(Collection< ? > c) {
		throw new UnsupportedOperationException();
	}

	public void clear() {
		throw new UnsupportedOperationException();
	}

	public Float get(int index) {
		return new Float(data[index]);
	}

	public float getFloat(int index) {
		return data[index];
	}

	public Float set(int index, Float element) {
		Float old = get(index);
		data[index] = element.floatValue();
		return old;
	}

	public void setFloat(int index, float element) {
		if (index > lastElement) {
			throw new ArrayIndexOutOfBoundsException();
		}
		data[index] = element;
	}
	
	public void add(int index, Float element) {
		setFloat(index, element.floatValue());		
	}

	public void addFloat(int index, float element) {
		setFloat(index, element);
	}

	public Float remove(int index) {
		throw new UnsupportedOperationException();
	}

	public int indexOf(Object o) {
		throw new UnsupportedOperationException();
	}

	public int lastIndexOf(Object o) {
		throw new UnsupportedOperationException();
	}

	public ListIterator<Float> listIterator() {
		throw new UnsupportedOperationException();
	}

	public ListIterator<Float> listIterator(int index) {
		throw new UnsupportedOperationException();
	}

	public List<Float> subList(int fromIndex, int toIndex) {
		throw new UnsupportedOperationException();
	}

	public float max() {
		float max = Float.MIN_VALUE;
		for (int i = 0; i < lastElement; i++) {
			if (data[i] > max) {
				max = data[i];
			}
		}
		return max;
	}

	public float min() {
		float min = Float.MAX_VALUE;
		for (int i = 0; i < lastElement; i++) {
			if (data[i] < min) {
				min = data[i];
			}
		}
		return min;
	}
	
	public double[] convertToPrimitiveDoubles() {
		double[] d = new double[size()];
		for (int i = 0; i < size(); i++) {
			d[i] = (double)data[i];
		}
		return d;
	}
	

	public Float[] convertToFloats() {
		Float[] f = new Float[size()];
		for (int i = 0; i < size(); i++) {
			f[i] = new Float(data[i]);
		}
		return f;		
	}

	public Double[] convertToDoubles() {
		Double[] d = new Double[size()];
		for (int i = 0; i < size(); i++) {
			d[i] = new Double(data[i]);
		}
		return d;		
	}
	
	public String toString() {
		String s = "";
		for (int i = 0; i < lastElement; i++) {
			s += data[i] + " ";
		}
		return s;
	}
}