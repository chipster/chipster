package fi.csc.microarray.databeans;

import java.util.List;

import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;

public class LinkUtils {

	/**
	 * <p>Retrieve a DataBean that is linked to given bean or to its ancestors. 
	 * Link type must be defined. If multiple ancestors have incoming
	 * link of the defined type, the closest one is chosen. If there are multiple
	 * incoming links, null is returned.</p>
	 * 
	 * <p>For example, by defining ANNOTATION link type we can retrieve
	 * the annotating dataset so that annotations are inherited.</p>
	 *  
	 */
	public static DataBean retrieveInherited(DataBean input, final Link linkType) {
		
		List<DataBean> inheritedBeans = input.traverseLinks(Link.derivationalTypes(), Traversal.DIRECT, new DataBeanSelector() {
			
			public boolean shouldSelect(DataBean bean) {				
				return bean.getLinkSources(linkType).size() > 0; // check if it is annotation
			}
			
			public boolean shouldTraverse(DataBean bean) {
				return true;
			}
		});
		
		if (inheritedBeans.size() > 0) {
			return inheritedBeans.get(0).getLinkSources(linkType).get(0);
		} else {
			return null;
		}
	}
}
