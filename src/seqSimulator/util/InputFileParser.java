package seqSimulator.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import seqSimulator.core.CoreSimulator;

public class InputFileParser {

	CoreSimulator m_coreSimulator;
	
	public CoreSimulator parseFile(File inputFile) throws Exception{
		
		parse(inputFile);
		
        if (m_coreSimulator != null)
            return m_coreSimulator;
        else {
            throw new Exception("CoreSimulator cannot be created from input file.");
        }
	}
	
	void parse(File inputFile) throws IOException{
	
		FileReader reader = new FileReader(inputFile);
		BufferedReader bufferedReader = new BufferedReader(reader);
		
	    String line = "";
	    while ((line = bufferedReader.readLine()) != null) {
	      System.out.println(line);
	    }
	    
	}
	
}
