package conifer.io.featurefactory;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

/** 
* @author Tingting Zhao (zhaott0416@gmail.com)
*/
public class JsonStringUtil {
	
	/**
	   * Convert a JSON string to pretty print version
	   * @param jsonString
	   * @return
	   */
	  public static String toPrettyFormat(String jsonString) 
	  {
	      JsonParser parser = new JsonParser();
	      JsonObject json = parser.parse(jsonString).getAsJsonObject();

	      Gson gson = new GsonBuilder().setPrettyPrinting().create();
	      String prettyJson = gson.toJson(json);

	      return prettyJson;
	  }
	  
	  
	  public static void writePrettyFormatToFile(Map<String, Set<String>> hashMap, String url){
		  // url saves the location and name of the file to be written 
		  // one example can be "/Users/crystal/Dropbox/backup/AminoAcidToCompressedCodons.txt"
		  String gson = new Gson().toJson(hashMap);
		  String prettyJson = toPrettyFormat(gson);
		  try(  PrintWriter out = new PrintWriter( url) ){
			    out.println( prettyJson );
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	  
		  
	  }
	  	 
	 
	  public static void writePrettyFormatToFile(String jsonString, String url){
		  String prettyJson = toPrettyFormat(jsonString);
		  try(  PrintWriter out = new PrintWriter( url) ){
			    out.println( prettyJson );
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	  		  
	  }
	  

}
