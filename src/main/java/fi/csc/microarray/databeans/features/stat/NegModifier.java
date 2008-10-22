package fi.csc.microarray.databeans.features.stat;


import fi.csc.microarray.databeans.features.CalculatingIterable;
import fi.csc.microarray.databeans.features.Modifier;

public class NegModifier extends CalcModifier implements Modifier {

	public NegModifier() {
		super(CalculatingIterable.CalcOperation.NEGATE);
	}

}
