package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.event.MouseEvent;

/**
 * Class that adds tooltip support to {@link GBrowserChartPanel}.
 * 
 * @author Aleksi Kallio
 *
 */
public class TooltipAugmentedChartPanel extends GBrowserChartPanel {

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
