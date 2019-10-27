package org.nupack.nudraw;

import java.util.Iterator;

import org.apache.log4j.Category;
import org.apache.log4j.Logger;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;

public class ConsoleRunner {

	private static final Category log = Logger.getLogger(ConsoleRunner.class);

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		JSAPResult config = parseArgs(args);
		
		DotParenParser parser = new DotParenParser(config.getString("inFile"));
		SequenceProcessor processor = new SequenceProcessor(parser.getParsedSequence());
		
		processor.runCalculationChain();
		
		log.debug("Number of strands: " + processor.getSequence().getNumStrands());
		log.debug("Number of bases: " + processor.getSequence().getNumBases());
		
		Iterator i = processor.getSequence().getBaseIterator();
		while(i.hasNext()) {
			log.debug(i.next());
		}
		
		
		return;
	}
	
	private static JSAPResult parseArgs(String[] args) {
		JSAP argParser = new JSAP();
		
		// set up argument parsing
		
		FlaggedOption inFileParam = new FlaggedOption("inFile")
									.setStringParser(JSAP.STRING_PARSER)
									.setRequired(true)
									.setShortFlag('i')
									.setLongFlag("in");
		inFileParam.setHelp("Input file in paren-dot format");
		
		FlaggedOption verbosityParam = new FlaggedOption("verbosity")
									.setStringParser(JSAP.INTEGER_PARSER)
									.setRequired(false)
									.setDefault("0")
									.setShortFlag('v')
									.setLongFlag("verbosity");
		verbosityParam.setHelp("Set verbosity level");
		
		Switch usageParam = new Switch("usage")
									.setShortFlag('u')
									.setLongFlag("usage");
		usageParam.setHelp("Display basic usage information");

		Switch helpParam = new Switch("help")
									.setShortFlag('h')
									.setLongFlag("help");
		usageParam.setHelp("Display help");
		
		try {
			argParser.registerParameter(inFileParam);
			argParser.registerParameter(verbosityParam);
			argParser.registerParameter(usageParam);
			argParser.registerParameter(helpParam);
		}
		catch(JSAPException e) {
			log.fatal("Couldn't add parameters to the arg parser: " + e.getMessage());
			System.exit(1);
		}
		
		JSAPResult config = argParser.parse(args);
		
		/*for(int i = 0; i < config.getInt("verbosity"); i++)
		{
			log.debug("got an entry");
		}*/
		
		// handle switches for help and usage
		if(config.getBoolean("usage")) {
			System.err.println("Usage: java "
								+ ConsoleRunner.class.getName()
								+ argParser.getUsage());
			System.exit(0);
		}
		
		if(config.getBoolean("help")) {
			System.err.print(argParser.getHelp());
			System.exit(0);
		}
		
		// if it some args were messed up, tell the user why
		if(! config.success()) {
			for(java.util.Iterator errors = config.getErrorMessageIterator(); errors.hasNext();) {
				System.err.println("Error: " + errors.next());
			}
			System.out.println(argParser.getUsage());
			
			System.exit(1);
		}
		
		// if we got to this point, nothing went wrong and no terminating options were found (help, usage, etc)
		
		return config;
	}

}
