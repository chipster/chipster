package fi.csc.microarray.client.dataimport.trimmer;

public class NormalStringReplace extends DataTrimmingOperation {

	private String oldString;
	private String newString;
	private boolean matchWholeString;

	public NormalStringReplace(String oldString, String newString, boolean matchWholeString, int columnIndex) {
		super(columnIndex);
		this.oldString = oldString;
		this.newString = newString;
		this.matchWholeString = matchWholeString;
	}
	
	@Override
	public String doTrimming(String stringToTrim) {
		if(matchWholeString){
			if(stringToTrim.equals(oldString)){
				return newString;
			} else {
				return stringToTrim;
			}
		} else {
			return stringToTrim.replace(oldString, newString);
		}
	}
	
	@Override
	public String toString(){
		return "NormalStringReplace [oldString="+oldString+", newString="+newString+", matchWholeString="+matchWholeString+", column="+getColumnIndex()+"]";
	}

}
