package fi.csc.chipster.web.adminweb.data;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

import org.hibernate.criterion.Criterion;
import org.hibernate.criterion.Restrictions;

import com.vaadin.data.hbnutil.filter.ContainerFilter;

public class DateContainerFilter extends ContainerFilter {

	private final Date searchDateStart;
	private final Date searchDateEnd;

	public DateContainerFilter(Object propertyId, String dateString) throws NumberFormatException {
		super(propertyId);

		Calendar startCal = new GregorianCalendar();
		Calendar endCal = new GregorianCalendar();

		String[] spaceSplit = dateString.split(" ");
		String date = spaceSplit[0];

		String[] hyphenSplit = date.split("-");

		String yearString = hyphenSplit[0];
		
		int year = -1;
		int month = -1;
		int dayOfMonth = -1;

		try {
			year = Integer.parseInt(yearString);
			
			//Month is 0-based
			startCal.set(year, 0, 1, 0, 0, 0);
			endCal.set(year, 11, 31, 23, 59, 59);

			if (hyphenSplit.length > 1) {
				String monthString = hyphenSplit[1];
				month = Integer.parseInt(monthString) - 1;//Month is 0-based
				
				startCal.set(year, month, 1);
				endCal.set(year, month, 31);

				if (hyphenSplit.length > 2) {
					String dayString = hyphenSplit[2];
					dayOfMonth = Integer.parseInt(dayString);
					
					startCal.set(year,  month, dayOfMonth);
					endCal.set(year, month, dayOfMonth);
				}
			}
		} catch (NumberFormatException e) {
			throw new NumberFormatException(
					"Search term parsing failed. Use following format: " +
					"2010 to search for year, " +
					"2010-03 to search for month or " +
					"2010-03-12 to search for specific date.");
		}
		

		this.searchDateStart = startCal.getTime();
		this.searchDateEnd = endCal.getTime();
	}

	@Override
	public Criterion getFieldCriterion(String fullPropertyName) {
		return Restrictions.between(fullPropertyName, searchDateStart, searchDateEnd);
	}

	public static String getToday() {
		Date now = new Date();
		String today = new SimpleDateFormat("yyyy-MM-dd").format(now);
		return today;
	}
}
