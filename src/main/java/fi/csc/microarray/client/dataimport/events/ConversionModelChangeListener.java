package fi.csc.microarray.client.dataimport.events;

import java.util.EventListener;

/**
 * Listener for changes in conversion model
 * 
 * @author mkoski
 *
 */
public interface ConversionModelChangeListener extends EventListener {

	public void headerChanged(HeaderChangedEvent e);
	public void footerChanged(FooterChangedEvent e);
	public void titleRowChanged(TitleRowChangedEvent e);
	public void delimiterChanged(DelimiterChangedEvent e);
	public void decimalSeparatorChanged(DecimalSeparatorChangedEvent e);
	public void columnTitlesChanged(ColumnTitlesChangedEvent e);
	public void inputFileChanged(InputFileChangedEvent e);
}
