package org.nupack.nudraw;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.log4j.Category;
import org.apache.log4j.Logger;

/**
 * Parses the nudraw input file format into the appropriate Base objects
 * 
 * @author mbp
 *
 */
public class DotParenParser {

	private static final char NICK_SYMBOL 			= '+';
	private static final char PAIR_LEFT_SYMBOL 		= '(';
	private static final char PAIR_RIGHT_SYMBOL 	= ')';
	private static final char UNPAIRED_SYMBOL 		= '.';
	private static final String SYMBOL_PATTERN 		= "^["
														+ NICK_SYMBOL
														+ PAIR_LEFT_SYMBOL
														+ PAIR_RIGHT_SYMBOL
														+ UNPAIRED_SYMBOL
														+ "]+$";
	private static final String COMMENT_PATTERN 	= "^\\s*[#%].*";
	
	private String filename = null;
	private Sequence sequence = null;
	private static final Category log = Logger.getLogger(DotParenParser.class);

	public DotParenParser(String filename) {
		this.filename = filename;
		this.sequence = new Sequence();
	}

	/* public interface to do something useful */
	public Sequence getParsedSequence() {
		parseInput();
		
		return sequence;
	}
	
	/**
	 * Read the input string and generate a sequence of strands & pairs
	 *
	 */
	private void parseInput() {
		String input = readInput();
		
		Base b = null;
		
		int baseIndex = 0;
		for (int charIndex = 0; charIndex < input.length(); charIndex++) {
			b = new Base();
			b.setIndex(baseIndex);
			
			switch (input.charAt(charIndex)) {
			// We don't add a symbol if it's a nick
			case NICK_SYMBOL:
				sequence.newStrand();
				continue;
			case UNPAIRED_SYMBOL:
				b.setType(BaseType.UNPAIRED);
				baseIndex++;
				break;
			case PAIR_LEFT_SYMBOL:
				b.setType(BaseType.PAIR_LEFT);
				baseIndex++;
				break;
			case PAIR_RIGHT_SYMBOL:
				b.setType(BaseType.PAIR_RIGHT);
				baseIndex++;
				break;
			default:
				throw new Error("Incorrect base symbol: "
						+ input.charAt(charIndex));
			}
			           
			sequence.addBase(b);
		}
	}


	/**
	 * Read the input file (ignoring comments, etc)
	 * @return Concatenated input
	 */
	private String readInput() {
		BufferedReader reader = null;
		String line = null;

		String allInput = "";

		try {
			reader = new BufferedReader(new FileReader(filename));

			while ((line = reader.readLine()) != null) {
				if (line.matches(COMMENT_PATTERN)) {
					continue;
				} else if (line.matches(SYMBOL_PATTERN)) {
					allInput += line;
				} else {
					log.info("Ignoring invalid input line: <" + line + ">");
				}
			}

		} catch (FileNotFoundException e) {
			log.fatal("Cannot open input file");
			System.exit(1);
		} catch (IOException e) {
			log.fatal("Error reading input file: " + e.getMessage());
			System.exit(1);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					log.fatal("Couldn't close the input stream");
					System.exit(1);
				}
			}
		}
		
		return allInput;
	}

}
