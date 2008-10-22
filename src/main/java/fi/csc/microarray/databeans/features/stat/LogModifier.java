package fi.csc.microarray.databeans.features.stat;


import fi.csc.microarray.databeans.features.CalculatingIterable;
import fi.csc.microarray.databeans.features.Modifier;

public class LogModifier extends CalcModifier implements Modifier {

	public LogModifier() {
		super(CalculatingIterable.CalcOperation.LOG_2);
	}

}
