package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SelectionText;

public class CnaRow implements SelectionText {
	
	public static class Sample {
		
		private Float flag;
		private Float logRatio;
		private String name;
		
		public Float getFlag() {
			return flag;
		}
		public void setFlag(Float flag) {
			this.flag = flag;
		}
		public Float getLogRatio() {
			return logRatio;
		}
		public void setLogRatio(Float logRatio) {
			this.logRatio = logRatio;
		}
		public String getName() {
			return name;
		}
		public void setName(String name) {
			this.name = name;
		}
		
		@Override 
		public String toString() {
			return "  " + getName() + 
					"\n   flag \t" + getFlag() + 
					"\n   segmented \t" + getLogRatio();
		}
	}
	
	private Region region;
	private LinkedList<Sample> samples;
	private Float lossFreg;
	private Float gainFreg;
	
	public Region getRegion() {
		return region;
	}
	public void setRegion(Region region) {
		this.region = region;
	}
	public LinkedList<Sample> getSamples() {
		return samples;
	}
	public void setSamples(LinkedList<Sample> samples) {
		this.samples = samples;
	}
	public Float getLossFreg() {
		return lossFreg;
	}
	public void setLossFreg(Float lossFreg) {
		this.lossFreg = lossFreg;
	}
	public Float getGainFreg() {
		return gainFreg;
	}
	public void setGainFreg(Float gainFreg) {
		this.gainFreg = gainFreg;
	}
	
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		
		builder.append("CNA feature info" + 
		"\n chromosome \t" + getRegion().start.chr + 
		"\n start \t" + getRegion().start.bp + 
		"\n end \t" + getRegion().end.bp + 
		"\n gain.freq \t" + getGainFreg() + 
		"\n loss.freq \t" + getLossFreg() + "\n");
				
		builder.append(" Samples\n");
		
		for (Sample sample : getSamples()) {			
			builder.append(sample);
			builder.append("\n");
		}
		
		return builder.toString();
	}
	
	@Override
	public String getText() {
		return this.toString();
	}
}
