package fi.csc.microarray.databeans.features;

import java.util.Iterator;

public class CalculatingIterable implements Iterable<Float> {

	public static enum CalcOperation {
		ADD,
		SUBTRACT,
		MULTIPLY,
		DIVIDE,
		LOG_2, 
		NEGATE;
	}

	private Iterable<Float> f1;
	private Iterable<Float> f2;
	private CalcOperation operation;

	public CalculatingIterable(Iterable<Float> f1, Iterable<Float> f2, CalcOperation operation) {
		this.f1 = f1;
		this.f2 = f2;
		this.operation = operation;
	}
	
	public static class CalculatingIterator implements Iterator<Float> {

		private CalcOperation operation;
		private Iterator<Float> f1;
		private Iterator<Float> f2;


		public CalculatingIterator(Iterator<Float> iterator, Iterator<Float> iterator2, CalcOperation operation) {
			this.f1 = iterator;
			this.f2 = iterator2;
			this.operation = operation;
		}

		public boolean hasNext() {
			return f1.hasNext() && (f2 != null ? f2.hasNext() : true);
		}

		public Float next() {
			return nextFloat();
		}

		public float nextFloat() {
			float r;
			switch (operation) {
			case ADD:
				r = f1.next() + f2.next();
				break;
			case SUBTRACT:
				r = f1.next() - f2.next();
				break;
			case MULTIPLY:
				r = f1.next() * f2.next();
				break;
			case DIVIDE:
				r = f1.next() / f2.next();
				break;
			case LOG_2:
				r = (float)(Math.log(f1.next()) / Math.log(2f)); // 2-based logarithm: log_2(x) = log_e(x) / log_e(2)
				break;
			case NEGATE:
				r = -f1.next(); // 2-based logarithm: log_2(x) = log_e(x) / log_e(2)
				break;
			default:
				throw new UnsupportedOperationException("unknown operation " + operation);
			}
			return r;
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
	
	public Iterator<Float> iterator() {
		Iterator<Float> iterator1 = f1.iterator();
		Iterator<Float> iterator2 = f2 != null ? f2.iterator() : null;
		return new CalculatingIterator(iterator1, iterator2, operation);
	}

}
