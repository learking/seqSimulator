package seqSimulator.app;

import java.io.File;
import java.util.Arrays;

import seqSimulator.core.CoreSimulator;
import seqSimulator.util.InputFileParser;

public class SeqSimulator {
	
	CoreSimulator m_coreSimulator;
	
    /**
     * parse command line arguments, and load file if specified
     *
     * @throws Exception *
     */
    public void parseArgs(String[] args) throws Exception {
        //int i = 0;
        File inputFile;
        
        try {
            //while (i < args.length) {
            	
                if (args.length == 1) {
                        inputFile = new File(args[0]);
                } else {
                        throw new Exception("Wrong argument");
                }
                
            //}
        }catch (Exception e) {
            e.printStackTrace();
            throw new Exception("Error parsing command line arguments: " + Arrays.toString(args) + "\nArguments ignored\n\n" + getUsage());
        }
        
        m_coreSimulator = new InputFileParser().parseFile(inputFile);
        
    }
	
    public static String getUsage() {
    	/*
        return "Usage: BeastMCMC [options] <Beast.xml>\n" +
                "where <Beast.xml> the name of a file specifying a Beast run\n" +
                "and the following options are allowed:\n" +
                "-resume : read state that was stored at the end of the last run from file and append log file\n" +
                "-overwrite : overwrite existing log files (if any). By default, existing files will not be overwritten.\n" +
                "-seed [<int>|random] : sets random number seed (default 127), or picks a random seed\n" +
                "-threads <int> : sets number of threads (default 1)\n" +
                "-prefix <name> : use name as prefix for all log files\n" +
                "-beastlib <path> : Colon separated list of directories. All jar files in the path are loaded. (default 'beastlib')";
 */
    	return "Usage: SeqSimulator [options] <>\n";
    } // getUsage
    
    public void run() throws Exception {
    	m_coreSimulator.run();
    }
    
    public static void main(String[] args) {
    	try{
    		SeqSimulator app = new SeqSimulator();
    		app.parseArgs(args);
    		app.run();
        } catch (Exception e) {
        	
        }
    	
    }
}