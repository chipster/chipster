package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedList;

public class CnaRow {
	
	public static class Sample {
		
		private float flag;
		private float logRatio;
		private String name;
		
		public float getFlag() {
			return flag;
		}
		public void setFlag(float flag) {
			this.flag = flag;
		}
		public float getLogRatio() {
			return logRatio;
		}
		public void setLogRatio(float logRatio) {
			this.logRatio = logRatio;
		}
		public String getName() {
			return name;
		}
		public void setName(String name) {
			this.name = name;
		}
	}
			
	private LinkedList<Sample> samples;
	private float lossFreg;
	private float gainFreg;
	
	public LinkedList<Sample> getSamples() {
		return samples;
	}
	public void setSamples(LinkedList<Sample> samples) {
		this.samples = samples;
	}
	public float getLossFreg() {
		return lossFreg;
	}
	public void setLossFreg(float lossFreg) {
		this.lossFreg = lossFreg;
	}
	public float getGainFreg() {
		return gainFreg;
	}
	public void setGainFreg(float gainFreg) {
		this.gainFreg = gainFreg;
	}
}
