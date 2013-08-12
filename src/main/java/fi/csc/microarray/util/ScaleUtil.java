package fi.csc.microarray.util;

import java.text.NumberFormat;

public class ScaleUtil {
	
	private static NumberFormat numberFormat = NumberFormat.getInstance();
	public static final int DEFAULT_STEP_COUNT = 4;
		
	public static String format(float f) {
		numberFormat.setMaximumFractionDigits(2);
		return numberFormat.format(f);
	}
	
	public static float[] generateScaleValues(float min, float max) {
		return generateScaleValues(min, max, DEFAULT_STEP_COUNT);
	}
	
	public static float[] generateScaleValues(float min, float max, int stepCount) {
		

		float[] scaleMinMax = getStepDistance(min, max, stepCount);
		
		//Scaling fails later if min equals max
		if (scaleMinMax[0] == scaleMinMax[1]) {
			scaleMinMax[1]++;
		}
		
		float scaleMin = scaleMinMax[0];
		float scaleMax = scaleMinMax[1];
		
		//If zero is in the scale, use it as second or third step of the scale.
		//This is implemented for 4 scale steps only
		if (scaleMin < 0 && scaleMax > 0 && stepCount == DEFAULT_STEP_COUNT) {
			
			if (Math.abs(min) > max) {
				//Zero on third step
				if (Math.abs(scaleMin) > scaleMax * 2) {
					scaleMax = Math.abs(scaleMin / 2);
				} else {
					scaleMin = -scaleMax * 2;
				}
			} else {
				//Zero on second step
				if (Math.abs(scaleMin) * 2 > scaleMax) {
					scaleMax = Math.abs(scaleMin * 2);
				} else {
					scaleMin = -scaleMax / 2;
				}
			}
		}

		float[] values = new float[stepCount];
		for (int i = 0; i < stepCount; i++) {
			values[i] = scaleMin + (scaleMax - scaleMin)
					/ (stepCount - 1) * i;
		}

		return values;
	}

	private static float[] getStepDistance(float minValue, float maxValue, int stepCount) {

		float preferredScale = Math.abs(maxValue - minValue) / (stepCount);
		float decimalFactor = calculateDecimalFactor(preferredScale);
		float scaleMax = (float) (Math.ceil(maxValue / decimalFactor))
				* decimalFactor;
		float scaleMin = (float) (Math.floor(minValue / decimalFactor))
				* decimalFactor;		

		return new float[] { scaleMin, scaleMax };
	}

	/**
	 * Calculates the number which is power of ten and just below preferredScale
	 * e.g. 1.3 -> 1, 23 -> 10, 543 -> 100
	 * 
	 * @param maxValue
	 * @return
	 */
	private static float calculateDecimalFactor(float preferredScale) {
		// To scale the preset SCALES to the size of the prefferredScale
		float decimalFactor = 1;
		// radix is 10
		if (preferredScale >= 10) {
			while (preferredScale >= decimalFactor * 10) {
				decimalFactor *= 10;
			}
		} else {
			while (preferredScale < decimalFactor / 10) {
				decimalFactor /= 10;
			}
		}
		return decimalFactor;
	}
}
