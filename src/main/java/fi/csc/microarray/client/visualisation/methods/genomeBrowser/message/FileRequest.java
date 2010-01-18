package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.TreeNode;

	public class FileRequest{
		
		
		

		public FileRequest(AreaRequest areaRequest, RowRegion rowRegion, TreeNode node, FsfStatus status) {
			super();
			this.rowRegion = rowRegion;
			this.node = node;
			this.status = status;
			this.areaRequest = areaRequest;
		}
		
		public AreaRequest areaRequest;
		public RowRegion rowRegion;
		public TreeNode node;
		
		public FsfStatus status;
	}