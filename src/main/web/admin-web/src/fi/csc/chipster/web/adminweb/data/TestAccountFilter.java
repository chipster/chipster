package fi.csc.chipster.web.adminweb.data;

import java.util.LinkedList;
import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.context.internal.ThreadLocalSessionContext;
import org.hibernate.criterion.Restrictions;

import com.vaadin.data.hbnutil.StringContainerFilter;

import fi.csc.chipster.web.adminweb.hbncontainer.HibernateUtil;

public class TestAccountFilter {

	private List<String> testAccounts;

	private List<String> getTestAccounts(Session session) {

		if (testAccounts == null) {
			@SuppressWarnings("unchecked")
			List<AccountEntry> results = session.createCriteria(AccountEntry.class)
			//FIXME
			//.add(Restrictions.eq("ignoreInStatistics", true))
			.list();

			testAccounts = new LinkedList<String>();

			for (AccountEntry entry : results) {
				testAccounts.add(entry.getUsername());
			}
		}
		return testAccounts;
	}

	public void addCriteriaForTestAccounts(Session session,
			boolean ignoreTestAccounts, Criteria criteria) {

		if (ignoreTestAccounts) {
			for (String account : getTestAccounts(session)) {
				criteria.add(Restrictions.not(Restrictions.eq(JobLogContainer.USERNAME, account)));
			}
		}
	}

	public String getHqlForTestAccounts(Session session, boolean ignoreTestAccounts) {

		String hql = "where not ";

		if (ignoreTestAccounts) {

			boolean isFirst = true;

			List<String> accounts = getTestAccounts(session);			

			for (String account : accounts) {

				//Check that we don't add anything dangerous to the query (only alphanumeric and not empty)
				if (account.matches("^[a-zA-Z0-9]+$")) {

					if (!isFirst) {
						hql += "and not ";
					}
					isFirst = false;

					hql += JobLogContainer.USERNAME + "='" + account + "' ";
				}
			}

			if (accounts.size() > 0) {

				return hql;
			}
		}

		return "";
	}
}