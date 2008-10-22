package fi.csc.microarray.client.dataimport.events;


public interface ConversionModelChangeSupport {

	public void addConversionChangeListener(ConversionModelChangeListener l);
	public void removeConversionChangeListener(ConversionModelChangeListener l);
	public void fireDecimalSeparatorChangeEvent(DecimalSeparatorChangedEvent event);
	public void fireDelimiterChangeEvent(DelimiterChangedEvent event);
	public void fireFooterChangeEvent(FooterChangedEvent event);
	public void fireHeaderChangeEvent(HeaderChangedEvent event);
	public void fireTitleRowChangeEvent(TitleRowChangedEvent event);
	public void fireColumnTitlesChangeEvent(ColumnTitlesChangedEvent event);
	public void fireInputFileChangeEvent(InputFileChangedEvent event);
	
}
