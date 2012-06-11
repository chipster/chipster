package fi.csc.microarray.manager.web.data;

import java.io.Serializable;

public class SearchFilter implements Serializable {
 
      private final String term;
      private final Object propertyId;
 
      public SearchFilter(Object propertyId, String searchTerm) {
            this.propertyId = propertyId;
            this.term = searchTerm;
      }

	public String getTerm() {
		return term;
	}

	public Object getPropertyId() {
		return propertyId;
	}
}