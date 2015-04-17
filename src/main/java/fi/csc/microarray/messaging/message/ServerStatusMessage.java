package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.SystemMonitorUtil;

/** 
 * @author klemela
 *
 */
public class ServerStatusMessage extends ChipsterMessage {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(ServerStatusMessage.class);		
	
	public static final String KEY_CPU_LOAD = "load";
	public static final String KEY_CPU_CORES = "cores";
	public static final String KEY_CPU_PERCENTS = "cpuPercents";
	public static final String KEY_MEM_USED = "memoryUsed";
	public static final String KEY_MEM_TOTAL = "memoryTotal";
	public static final String KEY_MEM_PERCENTS = "memoryPercents";
	public static final String KEY_DISK_USED = "diskUsed";
	public static final String KEY_DISK_TOTAL = "diskTotal";
	public static final String KEY_DISK_PERCENTS = "diskSpace";	
	public static final String KEY_SCHEDULED_JOBS = "scheduledJobs";
	public static final String KEY_RUNNING_JOBS = "runningJobs";
	public static final String KEY_HOST = "host";
	public static final String KEY_HOST_ID = "hostId";
	public static final String KEY_STATUS = "status";
	
	private double load;
	private int cores;
	private int cpuPercents;
	private long memUsed;
	private long memTotal;
	private int memPercents;
	private long diskUsed;
	private long diskTotal;
	private int diskPercents;
	private int scheduledJobs;
	private int runningJobs;
	private String host;
	private String hostId;
	private String status;
	
	public ServerStatusMessage() {
		// used by ActiveMq
	}
	
	public ServerStatusMessage(double load, int cores, int cpuPercents,
			long memUsed, long memTotal, int memPercents, long diskUsed,
			long diskTotal, int diskPercents) {
		
		this.load = load;
		this.cores = cores;
		this.cpuPercents = cpuPercents;
		this.memUsed = memUsed;
		this.memTotal = memTotal;
		this.memPercents = memPercents;
		this.diskUsed = diskUsed;
		this.diskTotal = diskTotal;
		this.diskPercents = diskPercents;
		
	}

	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		
		this.load = from.getDouble(KEY_CPU_LOAD);
		this.cores = from.getInt(KEY_CPU_CORES);
		this.cpuPercents = from.getInt(KEY_CPU_PERCENTS);
		this.memUsed = from.getLong(KEY_MEM_USED);
		this.memTotal = from.getLong(KEY_MEM_TOTAL);
		this.memPercents = from.getInt(KEY_MEM_PERCENTS);
		this.diskUsed = from.getLong(KEY_DISK_USED);		
		this.diskTotal = from.getLong(KEY_DISK_TOTAL);
		this.diskPercents = from.getInt(KEY_DISK_PERCENTS);
		this.scheduledJobs = from.getInt(KEY_SCHEDULED_JOBS);
		this.runningJobs = from.getInt(KEY_RUNNING_JOBS);
		this.host = from.getString(KEY_HOST);
		this.hostId = from.getString(KEY_HOST_ID);
		this.status = from.getString(KEY_STATUS);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
				
		mapMessage.setDouble(KEY_CPU_LOAD, this.load);
		mapMessage.setInt(KEY_CPU_CORES, this.cores);
		mapMessage.setInt(KEY_CPU_PERCENTS, this.cpuPercents);
		mapMessage.setLong(KEY_MEM_USED, this.memUsed);
		mapMessage.setLong(KEY_MEM_TOTAL,  this.memTotal);
		mapMessage.setInt(KEY_MEM_PERCENTS, this.memPercents);
		mapMessage.setLong(KEY_DISK_USED, this.diskUsed);
		mapMessage.setLong(KEY_DISK_TOTAL, this.diskTotal);			
		mapMessage.setInt(KEY_DISK_PERCENTS, this.diskPercents);
		mapMessage.setInt(KEY_SCHEDULED_JOBS, this.scheduledJobs);
		mapMessage.setInt(KEY_RUNNING_JOBS, this.runningJobs);
		mapMessage.setString(KEY_HOST, this.host);
		mapMessage.setString(KEY_HOST_ID, this.hostId);
		mapMessage.setString(KEY_STATUS, this.status);
	}
	
	public String toString() {
		String string = "";
		if (status != null) {
			string += status + "\n";			
		}
		string += "Jobs scheduled   \t" + scheduledJobs + "\n";
		string += "Jobs running     \t" + runningJobs + " \n";

		string += systemStatsToString();
		
		return string;
	}
	
	public String systemStatsToString() {
		String string = "";
		string += "Cpu %            \t" + cpuPercents +  "\t(load " + load +  ", cores " + cores +  ")\n";
		string += "Memory %         \t" + memPercents +  "\t(" + SystemMonitorUtil.bytesToGigas(memUsed) +  " / " + SystemMonitorUtil.bytesToGigas(memTotal) + " GB)\n";
		string += "Disk space %     \t" + diskPercents + "\t(" + SystemMonitorUtil.bytesToGigas(diskUsed) + " / " + SystemMonitorUtil.bytesToGigas(diskTotal) + " GB)\n";
		return string;
	}

	public String toStringLine() {
		String[] array = new String[] {				
				"" + scheduledJobs,
				"" + runningJobs,
				"" + cpuPercents,
				"" + load,
				"" + cores,
				"" + memPercents,
				SystemMonitorUtil.bytesToGigas(memUsed),
				SystemMonitorUtil.bytesToGigas(memTotal),
				"" + diskPercents,
				SystemMonitorUtil.bytesToGigas(diskUsed),
				SystemMonitorUtil.bytesToGigas(diskTotal),
				status == null ? "" : status
		};
		return StringUtils.rightPad(hostId, 40) + StringUtils.rightPad(host, 30) + Strings.rightPad(array, 10);
				
	}
	
	public static String getStringLineHeaders() {
		String[] array = new String[] {
				"SCHEDULED",
				"RUNNING",
				"CPU %",
				"LOAD",
				"CORES",
				"MEM %",
				"USED GB",
				"TOTAL GB",
				"DISK %",
				"USED GB",
				"TOTAL GB",
				"STATUS"
				};
		
		return StringUtils.rightPad("COMP ID", 40) + StringUtils.rightPad("HOST", 30) + Strings.rightPad(array, 10);		
	}

	public long getDiskPercents() {
		return diskPercents;
	}

	public void setDiskPercents(int diskPercents) {
		this.diskPercents = diskPercents;
	}

	public long getDiskTotal() {
		return diskTotal;
	}

	public void setDiskTotal(long diskTotal) {
		this.diskTotal = diskTotal;
	}

	public long getDiskUsed() {
		return diskUsed;
	}

	public void setDiskUsed(long diskUsed) {
		this.diskUsed = diskUsed;
	}

	public long getMemPercents() {
		return memPercents;
	}

	public void setMemPercents(int memPercents) {
		this.memPercents = memPercents;
	}

	public long getMemTotal() {
		return memTotal;
	}

	public void setMemTotal(long memTotal) {
		this.memTotal = memTotal;
	}

	public long getMemUsed() {
		return memUsed;
	}

	public void setMemUsed(long memUsed) {
		this.memUsed = memUsed;
	}

	public int getCpuPercents() {
		return cpuPercents;
	}

	public void setCpuPercents(int cpuPercents) {
		this.cpuPercents = cpuPercents;
	}

	public int getCores() {
		return cores;
	}

	public void setCores(int cores) {
		this.cores = cores;
	}

	public double getLoad() {
		return load;
	}

	public void setLoad(int load) {
		this.load = load;
	}

	public int getRunningJobs() {
		return runningJobs;
	}

	public void setRunningJobs(int runningJobs) {
		this.runningJobs = runningJobs;
	}

	public int getScheduledJobs() {
		return scheduledJobs;
	}

	public void setScheduledJobs(int scheduledJobs) {
		this.scheduledJobs = scheduledJobs;
	}

	public String getHost() {
		return this.host;
	}
	
	public void setHost(String host) {
		this.host = host;
	}
	
	public String getHostId() {
		return this.hostId;
	}

	public void setHostId(String id) {
		this.hostId = id;
	}

	public String getStatus() {
		return status;
	}

	public void setStatus(String status) {
		this.status = status;
	}
}
	

