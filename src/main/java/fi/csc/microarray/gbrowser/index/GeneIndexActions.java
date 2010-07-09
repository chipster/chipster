package fi.csc.microarray.gbrowser.index;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GeneIndexActions {

	private Connection conn;
	private Statement st;
	
	public GeneIndexActions(){
		
	}
	
    public void connect(){
        try {
			Class.forName("org.h2.Driver");
			conn = DriverManager.getConnection("jdbc:h2:~/test", "sa", "");
			/* 
			System.out.println("asdf");
			
			st.executeQuery("select * from testas");*/
			//Statement st = conn.createStatement();
			//st.executeUpdate("insert into testas values (null, 5)"); 
			
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("asdf "+e.getMessage());
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public void updateTable(List<RegionContent> indexList){
    	try {
			st = conn.createStatement();
			for (RegionContent id : indexList){
				st.executeUpdate("insert into index values (null,"+id.values.get(ColumnType.CHROMOSOME)+","+id.values.get(ColumnType.BP_START)+","+id.values.get(ColumnType.BP_END)+",'"+id.values.get(ColumnType.DESCRIPTION)+"' )");
    	}
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public GeneIndexDataType getLocation(String name){
    	try {
			st = conn.createStatement();
			ResultSet rs = st.executeQuery("select chromosome,bp_start,bp_end from index where name ='"+name+"'");
			rs.next();
			return new GeneIndexDataType(rs.getLong(1), rs.getLong(2), rs.getLong(3));
			
		} catch (SQLException e) {
			return null;
		}
    }
    
    public boolean checkIfNumber(String name){
    	try {
			Integer.parseInt(name);
			return true;
		} catch (NumberFormatException e) {
			try {
				Long.parseLong(name);
			} catch (NumberFormatException e1) {
				return false;
			}
			return false;
		}
    }
    
    public void closeConnection(){
    	try {
			conn.close();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
}
