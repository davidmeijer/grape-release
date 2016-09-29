package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.io.IOException;

import org.codehaus.jackson.map.ObjectMapper;

public class JsonUtil {

	
	public static void writeJson(String file, Object data) throws IOException {
		ObjectMapper mapper = new ObjectMapper();
		mapper.writeValue(new File(file), data);	
	}
}
