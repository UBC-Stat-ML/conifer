package conifer;


import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class RunFacility {
	
	public static String getCurrentDateString(){
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		//get current date time with Date()
		Date date = new Date();
		return dateFormat.format(date);
	}	
}
