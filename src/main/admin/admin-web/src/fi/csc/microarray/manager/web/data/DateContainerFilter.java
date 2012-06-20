package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

import org.hibernate.criterion.Criterion;
import org.hibernate.criterion.Restrictions;

import com.vaadin.data.hbnutil.ContainerFilter;

public class DateContainerFilter extends ContainerFilter {

	private final Date searchDateStart;
	private final Date searchDateEnd;

	public DateContainerFilter(Object propertyId, String dateString) {
		super(propertyId);

		Calendar startCal = new GregorianCalendar();
		Calendar endCal = new GregorianCalendar();

		String[] spaceSplit = dateString.split(" ");
		String fullDate = spaceSplit[0];
		//Don't care about the time 

		String[] hyphenSplit = fullDate.split("-");

		String yearString = hyphenSplit[0];
		
		int year = -1;
		int month = -1;
		int date = -1;

		try {
			year = Integer.parseInt(yearString);
			
			startCal.set(year, 1, 1);
			endCal.set(year, 12, 31);

			if (hyphenSplit.length > 1) {
				String monthString = hyphenSplit[1];
				month = Integer.parseInt(monthString);
				
				startCal.set(year, month, 1);
				endCal.set(year, month, 31);

				if (hyphenSplit.length > 2) {
					String dayString = hyphenSplit[2];
					date = Integer.parseInt(dayString);
					
					startCal.set(year,  month, date);
					endCal.set(year, month, date);
				}
			}
		} catch (NumberFormatException e) {
			//Parsing failed
		}

		this.searchDateStart = startCal.getTime();
		this.searchDateEnd = endCal.getTime();
	}

	@Override
	public Criterion getFieldCriterion(String fullPropertyName) {
		return Restrictions.between(fullPropertyName, searchDateStart, searchDateEnd);
	}
}
