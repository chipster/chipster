package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.event.MouseEvent;

import fi.csc.microarray.client.visualisation.NonScalableChartPanel;

public class TooltipEnabledChartPanel extends NonScalableChartPanel {

	public static interface TooltipRequestProcessor {
		public String tooltipRequest(MouseEvent mouseEvent);
	}
	
	private TooltipRequestProcessor tooltipRequestProcessor = null;

	public void addTooltipRequestProcessor(TooltipRequestProcessor tooltipRequestProcessor) {
		this.tooltipRequestProcessor = tooltipRequestProcessor;
	}
	
	public void removeTooltipRequestProcessor() {
		this.tooltipRequestProcessor = null;
	}
	
	@Override
	public String getToolTipText(MouseEvent me) {
		if (tooltipRequestProcessor != null) {
			return tooltipRequestProcessor.tooltipRequest(me);
			
		} else {
			return null;
		}
	}
}
