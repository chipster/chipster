package fi.csc.microarray.client.dataimport;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.ClosedByInterruptException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;

import sun.nio.ch.ChannelInputStream;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.events.ColumnTitlesChangedEvent;
import fi.csc.microarray.client.dataimport.events.ConversionModelChangeListener;
import fi.csc.microarray.client.dataimport.events.ConversionModelChangeSupport;
import fi.csc.microarray.client.dataimport.events.DecimalSeparatorChangedEvent;
import fi.csc.microarray.client.dataimport.events.DelimiterChangedEvent;
import fi.csc.microarray.client.dataimport.events.FooterChangedEvent;
import fi.csc.microarray.client.dataimport.events.HeaderChangedEvent;
import fi.csc.microarray.client.dataimport.events.InputFileChangedEvent;
import fi.csc.microarray.client.dataimport.events.TitleRowChangedEvent;

/**
 * Class to take care of the file reading/writing and chopping the data to
 * matrix using correct delimiter.
 * 
 * @author Mikko Koski, Petri Klemelä
 * 
 */
public class ConversionModel implements ConversionModelChangeSupport {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(ConversionModel.class);

	/**
	 * Column delimiter
	 */
	private Delimiter delim;

	/**
	 * Decimal separator
	 */
	private char decimalSeparator;

	/**
	 * Header end line number. If no header is set, the default value is
	 * negative value
	 */
	private int headerEnd = -1;

	/**
	 * Footer start line number. If no footer is set, the default value is an
	 * integer maximum value
	 */
	private int footerStart = Integer.MAX_VALUE;

	/**
	 * Column title line number. If column title line number is not set, the
	 * value is negative
	 */
	private int columnTitleLine = -1;

	/**
	 * Column titles. Includes also the first column title, which is practically
	 * an empty string. String starts with a number of column.
	 */
	private String[] columnTitles;
	
	
	/**
	 * Column titles without column number in the beginning to be used later as file names
	 */
	private String[] cleanColumnTitles;

	/**
	 * Listeners
	 */
	private List<ConversionModelChangeListener> conversionListeners;

	/**
	 * Input file
	 */
	private File inputFile;

	/**
	 * Count of lines on file. This is updated while file is analysed with
	 * <code>analyseFile</code>-method
	 */
	private int linesOnFile;

	/**
	 * Chopped data. Notice that the first column is a line number column
	 */
	private Object[][] choppedDataMatrix;

	/**
	 * Import screen
	 */
	private ImportScreenModel screen;

	/**
	 * Row limit. If the limit is not set, the default value is an integer
	 * maximum values.
	 */
	private int rowLimit = Integer.MAX_VALUE;

	/**
	 * Column limit. If the limit is not set, the default value is an integer
	 * maximum values. The first column, row number column, is also inclueded to
	 * this limit
	 */
	private int columnLimit = Integer.MAX_VALUE;

	/**
	 * The biggest number of columns in the row. The column count per line may
	 * vary because of the headers and footers
	 */
	private int maxColumnsInRow = 0;

	private FileAnalyser analyser;

	private BufferedReader buf;

	private List<File> outputFiles;

	/**
	 * The line of the file that reader is reading
	 */
	private int lineNumber;

	/**
	 * Store the count of lines read from file. The header and footer may be
	 * skipped so that's why the <code>lineNumber</code> does not necesserly
	 * represent the count of lines really chopped.
	 */
	private int linesRead;

	private ClientApplication application = Session.getSession().getApplication();

	/**
	 * Creates a conversion model without any limits
	 * 
	 * @param inputFile
	 *            File to be imported
	 */
	public ConversionModel(ImportScreenModel screen) {
		this(screen, Integer.MAX_VALUE, Integer.MAX_VALUE);

	}

	/**
	 * Creates a conversation model with limits
	 * 
	 * <strong>NOTICE!</strong> The column limit includes also the first column
	 * which is the row number column. So if you want to get 4 data columns the
	 * limit must be 5.
	 * 
	 * @param inputFile
	 *            File to be imported
	 * @param rowLimit
	 *            number of rows
	 * @param columnLimit
	 *            number of columns including the row number column
	 */
	public ConversionModel(ImportScreenModel screen, int rowLimit, int columnLimit) {
		this.rowLimit = rowLimit;
		this.columnLimit = columnLimit;
		this.screen = screen;
		this.headerEnd = -1;
		this.footerStart = Integer.MAX_VALUE;
		this.conversionListeners = new Vector<ConversionModelChangeListener>();
	}

	/**
	 * Gets decimal separator
	 * 
	 * @return decimal separator
	 */
	public char getDecimalSeparator() {
		return decimalSeparator;
	}

	/**
	 * Gets delimiter
	 * 
	 * @return delimiter
	 */
	public Delimiter getDelim() {
		return delim;
	}

	/**
	 * Sets delimiter and notifies listeners.
	 * 
	 * @param delim
	 *            new delimiter
	 */
	public void setDelim(Delimiter delim) {
		this.delim = delim;
		this.fireDelimiterChangeEvent(new DelimiterChangedEvent(this, delim));
	}

	/**
	 * Sets decimal sepataror and notifies listeners.
	 * 
	 * @param decimalSeparator
	 *            new decimal separator
	 */
	public void setDecimalSeparator(char decimalSeparator) {
		this.decimalSeparator = decimalSeparator;
		this.fireDecimalSeparatorChangeEvent(new DecimalSeparatorChangedEvent(this, decimalSeparator));
	}

	/**
	 * Gets line number where the header ends.
	 * 
	 * Default value if footer is not set: <strong>-1</strong>
	 * 
	 * @return header end line number or negative number if header is not set
	 */
	public int getHeaderEnd() {
		return headerEnd;
	}

	/**
	 * Gets line number where the footer starts.
	 * 
	 * Default value if footer is not set: <strong>Integer.MAX_VALUE</strong>
	 * 
	 * @return footer start line number an integer maximum value if footer is
	 *         not set
	 */
	public int getFooterStart() {
		return footerStart;
	}

	/**
	 * Gets default column titles which are increased number starting from one.
	 * This is used when user is not selected title row from the table.
	 * 
	 * @return default column titles (1,2,3...)
	 */
	private String[] getDefaultColumnTitles() {
		String[] columnTitles = new String[getLimitedColumnCount()];

		for (int i = 0; i < columnTitles.length; i++) {
			if (i == 0) {
				// Row number column
				// NOTICE! The whitespace is important. With out it, the columns
				// are displayed extremely small sized.
				columnTitles[0] = " ";
			} else {
				columnTitles[i] = "" + i;
			}
		}
		return columnTitles;
	}

	/**
	 * Analyses file, reads data and chops it at the same time.
	 * 
	 * This is the most important method of the whole data import process.
	 * 
	 * The data importing is started by first analysing the file to be imported.
	 * This is done with the help of <code>analyseFile</code> method, which
	 * reads the line count of the file and tries to guess the delimiter of the
	 * file.
	 * 
	 * The file is read line by line from the start to the end or just a limited
	 * amount of rows and columns. When a single line is read it is splitted and
	 * the row number column is added. Finally the chopped data is returned and
	 * saved to <code>choppedDataMatrix</code> field.
	 * 
	 * The method uses also <code>ProgressInformator</code> which is an
	 * interface to a informator class to inform the user of the process state.
	 * 
	 * @param ignoreHeadersAndFooters
	 *            ignores header and footer lines. This is used usually on the
	 *            second step of import tool
	 * @param informator
	 *            progress informator
	 * @throws IOException
	 * 
	 */
	public Object[][] chopData(boolean ignoreHeadersAndFooters, ProgressInformator informator) throws IOException {

		// Cannot read file and chop data if no input file is set
		if (inputFile == null) {
			throw new IllegalStateException("File can not be read if no input file is set");
		}

		informator.setIndeterminate(true);
		informator.setMessage("Analysing file");
		analyseInputFile(); // Throws IOException if interrupted

		// Sets minimum and maximum to informator
		informator.setMinimunValue(0);
		informator.setMaximumValue(this.getLimitedRowCount());
		informator.setIndeterminate(false);
		informator.setMessage("Parsing data and updating table");

		// Initialize variable which keeps track of the maximum column number
		// per row
		this.maxColumnsInRow = 0;

		// Gets start time (debugging)
		long chopStarted = System.currentTimeMillis();

		// LinkedList to store the splitted lines.
		List<Object[]> dataMatrix = new LinkedList<Object[]>();

		// Gets lines from file and split them on delimeter
		String[] splittedLine;
		Object[] dataLine;
		StringBuffer token;
		Double tokenAsDouble;

		informator.setMessage("Reading and parsing data from file");

		boolean isLimited = false;

		// Initialize buffered reader
		initializeReader();
		try {
			while ((splittedLine = nextLine(ignoreHeadersAndFooters, false)) != null && !Thread.interrupted()) {
				isLimited = false;

				// Sets the lenght of the line. The lenght is depending on
				// columnLimit
				int dataLineLenghtLimit;
				if (splittedLine.length >= columnLimit) {
					// Limit the lenght
					dataLineLenghtLimit = columnLimit;
					isLimited = true;
				} else {
					// No need to limit the lenght. Plus 1 for the row number
					// column
					dataLineLenghtLimit = splittedLine.length + 1;
				}

				// Count the actual column count
				if (splittedLine.length + 1 > this.maxColumnsInRow) {
					this.maxColumnsInRow = splittedLine.length + 1;
				}

				// Read the splitted line and store the information to dataLine
				// array
				dataLine = new Object[dataLineLenghtLimit];

				// Row number
				dataLine[0] = lineNumber;

				// Start from 1 because of the first column (row number)
				// Respects the column limit
				for (int j = 1; j < dataLineLenghtLimit; j++) {

					token = new StringBuffer(splittedLine[j - 1]);

					// Count the maximum column count in the file.
					// Plus one because we want to store the count, not the
					// index
					if (j + 1 > this.maxColumnsInRow) {
						this.maxColumnsInRow = j + 1;

						// If the maxColumnsInRow is changed and titles are the
						// default titles,
						// the change of column count affect also to column
						// titles.
						if (columnTitles == null && this.maxColumnsInRow == getLimitedColumnCount()) {
							this.fireColumnTitlesChangeEvent(new ColumnTitlesChangedEvent(this, getDefaultColumnTitles()));
						}
					}

					// To have visual sign of limitation
					if (isLimited && j == dataLineLenghtLimit - 1) {
						dataLine[j] = "...";
					} else {

						// Saves token to dataMatrix as Double or String
						try {
							tokenAsDouble = new Double(token.toString());
							dataLine[j] = tokenAsDouble;
						} catch (NumberFormatException nfe2) {
							dataLine[j] = token.toString();
						}
					}
				}

				// Add the line to the matrix
				dataMatrix.add(dataLine);

				// Inform the informator is needed
				if (lineNumber % getInformationInterval() == 0) {
					informator.setValue(lineNumber);
				}
			}
		} finally {
			// The buf.readline method may throw ClosedByInterruptException
			buf.close();
		}

		long chopTime = System.currentTimeMillis() - chopStarted;

		logger.debug("Chop time: " + chopTime);

		// Save the data to memory
		this.choppedDataMatrix = dataMatrix.toArray(new Object[0][]);

		return this.choppedDataMatrix;
	}

	/**
	 * Creates analyser if needed and analyses the file if it is not already
	 * analysed. Sets delimeter and linesOnFile after analysing is done
	 * 
	 * @throws IOException
	 */
	private void analyseInputFile() throws IOException {
		if (analyser == null) {
			// Creates analyser if it is not already created
			analyser = new FileAnalyser(inputFile);
		}

		if (analyser.isSameFileAs(inputFile)) {
			if (!analyser.isAnalysed()) {
				analyser.startAnalysing();
				this.setDelim(analyser.getDelimiter());
				this.setLinesOnFile(analyser.getLinesOnFile());
			}
		} else {
			analyser.setInputFile(inputFile);
			analyser.startAnalysing();
			this.setDelim(analyser.getDelimiter());
			this.setLinesOnFile(analyser.getLinesOnFile());
		}
	}

	private void setLinesOnFile(int linesOnFile) {
		this.linesOnFile = linesOnFile;
	}

	/**
	 * Gets information interval. The interval value is row count / 100 or 10,
	 * if the result of division is less than 10
	 * 
	 * @return information interval
	 */
	private int getInformationInterval() {
		int interval = getLimitedRowCount() / 100;
		if (interval < 10) {
			return 10;
		} else {
			return interval;
		}
	}

	private void initializeReader() throws FileNotFoundException {
		if (inputFile == null) {
			throw new IllegalStateException("Can't initialize reader while input file is null");
		}

		buf = new BufferedReader(new InputStreamReader(new ChannelInputStream((new RandomAccessFile(inputFile, "r")).getChannel())));
		lineNumber = 0;
		linesRead = 0;
	}

	/**
	 * The columns per row may vary because of the headers and footers. This
	 * method return the maximum number of columns in the row which contains the
	 * biggest count of columns in the whole table.
	 * 
	 * @return
	 */
	private int getMaximumColumnsInRow() {
		return maxColumnsInRow;
	}

	/**
	 * Gets the number of columns. If column limit is bigger than the actual
	 * column count, the limit is returned. Otherwise the actual column count is
	 * returned
	 * 
	 * @return columns on the matrix or column limit if the actual column count
	 *         is bigger than the limit.
	 */
	public int getLimitedColumnCount() {
		if (getMaximumColumnsInRow() > columnLimit) {
			return columnLimit;
		} else {
			return getMaximumColumnsInRow();
		}
	}

	/**
	 * Returns the total column count and ignores column limits
	 * 
	 * @return total column count
	 */
	public int getUnlimitedColumnCount() {
		return getMaximumColumnsInRow();
	}

	/**
	 * Returns the total row count (lines on file)
	 * 
	 * @return total row count (lines on file)
	 */
	public int getUnlimitedRowCount() {
		return linesOnFile;
	}

	/**
	 * Gets number of rows. If row limit is bigger than actual row count on
	 * file, the limit is returned
	 * 
	 * @return row on the matrix or row limit if the actual row count is bigger
	 *         than the limit.
	 */
	public int getLimitedRowCount() {
		if (linesOnFile > rowLimit) {
			return rowLimit;
		} else {
			return linesOnFile;
		}
	}

	/**
	 * Writes data to file. The file writing is done by reading the file line by
	 * line and then using conversion models and column managers setting the
	 * file is chopped and then written to file.
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public void writeToFile(ProgressInformator informator) throws FileNotFoundException, IOException {

		// Creates output files
		createOutputFiles(screen.getColumnTypeManager().getChipCount());

		// Initialises file reader
		initializeReader();

		// Creates writer for each chip (file)
		List<PrintWriter> outputWriters = new ArrayList<PrintWriter>();
		for (File outputFile : getOutputFiles()) {
			FileOutputStream outStream = new FileOutputStream(outputFile);
			outputWriters.add(new PrintWriter(new BufferedWriter(new OutputStreamWriter(outStream))));
		}

		// Initialize informator
		informator.setMinimunValue(0);
		informator.setMaximumValue(linesOnFile);

		String OUTPUT_DELIM = "\t";

		try {

			// Gets columns
			List<DataColumn> columns = screen.getColumnTypeManager().getColumns();

			boolean[] isFirst = new boolean[screen.getColumnTypeManager().getChipCount()];
			for (int i = 0; i < isFirst.length; i++) {
				isFirst[i] = true;
			}

			// Iterate through all columns and write column headers
			for (DataColumn column : columns) {
				ColumnType type = column.getColumnType();

				// Skip null, unused and row number columns
				if (type == null || type.equals(ColumnType.UNUSED_LABEL) || type.equals(ColumnType.ROW_NUMBER)) {
					continue;
				}

				if (type.equals(ColumnType.ANNOTATION_LABEL) || type.equals(ColumnType.IDENTIFIER_LABEL)) {
					for (int chip = 0; chip < outputWriters.size(); chip++) {
						// Do NOT print delimiter if it is first on the line
						if (!isFirst[chip]) {
							// Print delimeter
							outputWriters.get(chip).print(OUTPUT_DELIM);
						} else {
							// No delimeter before first value
							isFirst[chip] = false;
						}

						// Print column identifier
						outputWriters.get(chip).print(type.getIdentifier());
					}
				} else {
					// Do NOT print delimiter if it is first on the line
					if (!isFirst[column.getChipNumber() - 1]) {
						// Print delimeter
						outputWriters.get(column.getChipNumber() - 1).print(OUTPUT_DELIM);
					} else {
						// No delimeter before first value
						isFirst[column.getChipNumber() - 1] = false;
					}

					// Print column identifier
					outputWriters.get(column.getChipNumber() - 1).print(type.getIdentifier());
				}
			}

			// Finally, add new line to each file
			for (PrintWriter output : outputWriters) {
				output.print("\n");
			}

			/*
			 * Column identifiers are now written to each file, so let's write
			 * the actual data
			 */

			// Read the file line by line
			String[] splittedLine;
			while ((splittedLine = nextLine(true, true)) != null) {

				// Column number of the current line
				int columnNumber = 0;

				// Initializes isFirst variables to be true
				for (int i = 0; i < isFirst.length; i++) {
					isFirst[i] = true;
				}

				// Iterate through all columns
				for (DataColumn column : columns) {

					// Get column type
					ColumnType type = column.getColumnType();

					// Skip unused columns
					if (type.equals(ColumnType.UNUSED_LABEL)) {
						columnNumber++;
						continue;
					}

					// Skip row number column
					else if (type.equals(ColumnType.ROW_NUMBER)) {
						// The row number does not exists on the file, so no
						// need to increase columnNumber
						continue;
					}

					if (type.equals(ColumnType.ANNOTATION_LABEL) || type.equals(ColumnType.IDENTIFIER_LABEL)) {
						for (int chip = 0; chip < outputWriters.size(); chip++) {
							writeSplittedLineToChipFile(chip, columnNumber, isFirst, splittedLine, outputWriters.get(chip));
						}
					} else {
						writeSplittedLineToChipFile(column.getChipNumber() - 1, columnNumber, isFirst, splittedLine, outputWriters.get(column.getChipNumber() - 1));
					}

					columnNumber++;
				}

				// Add new line
				for (PrintWriter output : outputWriters) {
					output.print("\n");
				}

				if (lineNumber % getInformationInterval() == 0) {
					informator.setValue(lineNumber);
				}
			}
		} catch (ClosedByInterruptException cbie) {
			// Do nothing

		} catch (Exception e) {
			application.reportException(e);
		} finally {
			for (PrintWriter output : outputWriters) {
				if (output != null) {
					output.close();
				}
			}
			buf.close();
			choppedDataMatrix = null;
			informator.setValue(0);
		}
	}

	/**
	 * Sets the header end row
	 * 
	 * @param row
	 *            header end row
	 */
	public void setHeaderEndsRow(int row) {
		this.headerEnd = row;
		this.fireHeaderChangeEvent(new HeaderChangedEvent(this, row));
	}

	/**
	 * Sets the footer begin row
	 * 
	 * @param row
	 *            footer begin row
	 */
	public void setFooterBeginsRow(int row) {
		this.footerStart = row;
		this.fireFooterChangeEvent(new FooterChangedEvent(this, row));
	}

	public boolean hasHeader() {
		return this.getHeaderEnd() >= 0;
	}

	public boolean hasFooter() {
		return this.getFooterStart() != Integer.MAX_VALUE;
	}

	/**
	 * Sets the line which contains cell headers. After the line is set this
	 * method also read the line data and saves it so the column strings can be
	 * get with <code>getColumnTitles</code> method.
	 * 
	 * The first string of column titles is always set to white space because it
	 * is the row number column.
	 * 
	 * NOTICE! The columnTitleLine is the column line <strong>on the file</code>
	 * not on the table. If header and footer are not visible on the table the
	 * 23th row may be totally different that the 23th line on the table
	 * 
	 * The count of the titles is also limited after column limit
	 * 
	 * @param row
	 *            line which contains cell headers
	 */
	public void setColumnTitleLine(int row) {
		this.columnTitleLine = row;

		// Set titles
		if (row < 0) {
			this.columnTitles = null;
			logger.debug("Column titles set to null");
		} else {
			String[] columns = new String[getLimitedColumnCount()];
			for (int i = 1; i < columns.length && i < columnLimit; i++) {
				// The first column of choppedDataMatrix is ignored because it
				// is the
				// row number column
				if (i < choppedDataMatrix[row].length) {
					columns[i - 1] = choppedDataMatrix[row][i].toString();
				} else {
					columns[i - 1] = " ";
				}
			}
			this.setColumnTitles(columns);
			logger.debug("Column titles set. Column count is " + columnTitles.length);
		}

		this.fireTitleRowChangeEvent(new TitleRowChangedEvent(this, row));
	}

	/**
	 * Sets the column titles. Adds the empty title to the first column title
	 * and respects the column limit. This method does not notify listeners
	 * 
	 * @param columns
	 */
	private void setColumnTitles(String[] columns) {

		// Do not change the titles if there is nothing to change
		if (columns.equals(this.columnTitles)) {
			return;
		}

		// Create array for column titles. The lenght of array is
		// limited by columnLimit. The plus one comes from the first column.
		// columnLimit includes the first column.
		int columnTitlesLength;
		if (columns.length + 1 > columnLimit) {
			columnTitlesLength = columnLimit;
		} else {
			columnTitlesLength = columns.length + 1;
		}
		this.columnTitles = new String[columnTitlesLength];
		this.cleanColumnTitles = new String[columnTitlesLength];
		
		this.cleanColumnTitles[0] = " ";
		for (int i = 1; i < cleanColumnTitles.length; i++) {
			// Add column number to the beginning of the column name
			this.cleanColumnTitles[i] = columns[i - 1];
		}
				
		for (int i = 1; i < columnTitles.length; i++) {
			// Add column number to the beginning of the column name
			this.columnTitles[i] = "" + i + " - " + cleanColumnTitles[i];
		}

		fireColumnTitlesChangeEvent(new ColumnTitlesChangedEvent(this, this.columnTitles));
	}

	/**
	 * Sets row and column limits to the conversation model.
	 * 
	 * The column count may change if column limit is changed, so the column
	 * titles are also updated.
	 * 
	 * <strong>NOTICE!</strong> The column limit includes also the first column
	 * which is the row number column. So if you want to get 4 data columns the
	 * limit must be 5.
	 * 
	 * @param rowLimit
	 *            how many rows will be read from file
	 * @param columnLimit
	 *            how many columns will be read from file
	 */
	public void setLimits(int rowLimit, int columnLimit) {

		boolean columnChanged;
		if (this.columnLimit == columnLimit) {
			columnChanged = false;
		} else {
			columnChanged = true;
		}

		this.columnLimit = columnLimit;
		this.rowLimit = rowLimit;

		// If the column limit changes and the titles are not set, then the
		// limit
		// change affects to the columns too. Otherwise the columns are changed
		// after they have been read from the file in the nextLine() method
		if (columnChanged && columnTitles == null) {
			this.fireColumnTitlesChangeEvent(new ColumnTitlesChangedEvent(this, getDefaultColumnTitles()));
		}
	}

	/**
	 * Gets the column titles. If column titles is not set, the default column
	 * titles are returned. Default column titles are 1,2,3...
	 * 
	 * @return
	 */
	public String[] getColumnTitles() {
		if (this.columnTitles != null) {
			logger.debug("Got the predefined column titles. Column count is " + this.getLimitedColumnCount());
			return columnTitles;
		} else {
			logger.debug("Got default column titles. Column count is " + this.getLimitedColumnCount());
			return getDefaultColumnTitles();
		}
	}
	
	
	/**
	 * @return original column titles, null if not selected
	 */
	public String[] getCleanColumnTitles() {		
		return cleanColumnTitles;		
	}
	
	public String getCleanColumnTitle(int index) {
		if(cleanColumnTitles != null && cleanColumnTitles.length >= index) {
			return cleanColumnTitles[index];
		}
		return null;
	}
	

	/**
	 * Sets input file
	 * 
	 * @param inputFile
	 *            input file
	 */
	public void setInputFile(File inputFile) {
		this.inputFile = inputFile;

		fireInputFileChangeEvent(new InputFileChangedEvent(this, inputFile));
	}

	/**
	 * Retuns if the column title line is set
	 * 
	 */
	public boolean hasColumnTitles() {
		return !(columnTitleLine < 0);
	}

	/**
	 * Reads and returns next line of the file. The buffered reader must be
	 * initialized before using this method.
	 * 
	 * This method also stores the columnTitles if it reads the line which is
	 * marked to be column title line.
	 * 
	 * @param ignoreHeadersAndFooters
	 *            do not return header, footer or title row
	 * @param ignoreLimits
	 *            ignore row and column limits
	 * @return splitted line or <code>null</code> if end of file or limits
	 *         reached
	 * 
	 */
	private String[] nextLine(boolean ignoreHeadersAndFooters, boolean ignoreLimits) throws IOException, ClosedByInterruptException {
		String nextLine = null;
		while ((nextLine = buf.readLine()) != null && !Thread.interrupted()) {

			// Empty line
			if (nextLine.equals("")) {
				// Ignores empty lines and acts like they didn't even exists
				// lineNumber is not increased
				continue;
			}

			// Row limit
			else if (!ignoreLimits && linesRead >= rowLimit) {
				return null;
			}

			// Title row
			else if (lineNumber == columnTitleLine) {
				// Found the column line
				this.setColumnTitles(nextLine.split(delim.toString()));
				logger.debug("Column title line found. Count of columns: " + (nextLine.split(delim.toString())).length);
				lineNumber++;
				if (ignoreHeadersAndFooters) {
					continue;
				} else {
					linesRead++;
					return nextLine.split(delim.toString());
				}
			}

			// Header
			else if (ignoreHeadersAndFooters && (headerEnd > -1 && lineNumber <= headerEnd)) {
				lineNumber++;
				continue;
			}

			// Footer
			else if (ignoreHeadersAndFooters && (footerStart > -1 && lineNumber >= footerStart)) {
				lineNumber++;
				continue;
			}

			// Data
			else {
				lineNumber++;
				linesRead++;
				return nextLine.split(delim.toString());
			}
		}
		return null;
	}

	public void addConversionChangeListener(ConversionModelChangeListener l) {
		this.conversionListeners.add(l);
	}

	public void removeConversionChangeListener(ConversionModelChangeListener l) {
		this.conversionListeners.remove(l);
	}

	public void fireDecimalSeparatorChangeEvent(DecimalSeparatorChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.decimalSeparatorChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireDelimiterChangeEvent(DelimiterChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.delimiterChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireFooterChangeEvent(FooterChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.footerChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireHeaderChangeEvent(HeaderChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.headerChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireTitleRowChangeEvent(TitleRowChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.titleRowChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireColumnTitlesChangeEvent(ColumnTitlesChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.columnTitlesChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	public void fireInputFileChangeEvent(InputFileChangedEvent event) {
		for (ConversionModelChangeListener listener : conversionListeners) {
			listener.inputFileChanged(event);
		}
		logger.debug("Fired event " + event);
	}

	/**
	 * Creates output files. The output file's filename is based on input
	 * filename and chip number.
	 * 
	 * @param chipCount
	 *            count of chips in this file. This is also the count of output
	 *            file's since every chip is written on its own file
	 * @throws IOException
	 */
	public void createOutputFiles(int chipCount) throws IOException {
		if (outputFiles != null) {
			outputFiles.clear();
		} else {
			outputFiles = new ArrayList<File>();
		}

		//remove last file extension with a help of mystical regexp
		String originalFileName = inputFile.getName().replaceAll ("\\.[^.]*$", "");
		
		if (chipCount > 1) {
			
			List<String> columns = screen.getColumnTypeManager().getOriginalChipNames();
			
			for (int i = 0; i < columns.size(); i++) {
				String columnName = columns.get(i);
				if(columnName == null){
					columnName = originalFileName + "_" + (i + 1);
				}
				String fileName = columnName + " (";
				String suffix = ").tsv";
				outputFiles.add(ImportUtils.createTempFile(fileName, suffix));
			}
		} else {			
			String suffix = ").tsv";
			outputFiles.add(ImportUtils.createTempFile(originalFileName + " (", suffix));
		}
	}

	/**
	 * Gets the output files. The output files must be created by method
	 * <code>createOutputFiles</code>. If the files are not created before
	 * this method is called, returns <code>null</code>
	 * 
	 * @return output files
	 */
	public List<File> getOutputFiles() {
		return outputFiles;
	}

	private void writeSplittedLineToChipFile(int chip, int columnNumber, boolean[] isFirst, String[] splittedLine, PrintWriter outputWriter) {
		String OUTPUT_DELIM = "\t";

		if (isFirst[chip]) {
			// Do not print delimeter to the first column of the row
			isFirst[chip] = false;
		} else {
			outputWriter.print(OUTPUT_DELIM);
		}

		// Get the data to be written to file
		String dataToWrite;
		if (columnNumber >= splittedLine.length) {
			// Not enough column on this line.
			dataToWrite = "";
		} else {
			dataToWrite = splittedLine[columnNumber].toString();
		}

		// Do flag replacements ( +1 because there is no row number row)
		dataToWrite = screen.getFlagTrimmer().doTrimming(dataToWrite, columnNumber + 1);

		// Do data trimming replacements
		dataToWrite = screen.getDataTrimmer().doTrimming(dataToWrite, columnNumber + 1);

		// Get the correct outputWriter and print the data to file
		outputWriter.print(dataToWrite);
	}

	public String getInputFileName() {
		if (inputFile != null) {
			return inputFile.getName();
		} else {
			return null;
		}
	}
}
