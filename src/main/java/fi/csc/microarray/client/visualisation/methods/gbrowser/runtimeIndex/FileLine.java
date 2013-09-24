package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

public class FileLine {

	public String format(String title, String[] columns) {
		StringBuilder builder = new StringBuilder();
		builder.append(title);
		builder.append("\n");

		for (int i = 0; i < columns.length; i += 2) {
			String name = columns[i];
			String value = columns[i+1];

			//Unfortunately also null is converted to string sometimes
			if (!"null".equals(value)) {
				//builder.append("   ");
				builder.append(name);
				builder.append("\t");
				builder.append(value);
				builder.append("\n");
			}
		}

		return builder.toString();
	}
}
