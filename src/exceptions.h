#include <string>
#include <iostream>
#include <sstream>

/**
 * The basic class for exceptions in %Karlsberg.
 * Output to ostreams is enabled by the operator<<(std::ostream&,KbException&).
 * @author Bjoern Rabenstein
 */
class GeneralException 
{
public:
    /**
     * Constructor for a basic exception
     * @param errorCode the error code
     * @param message a descriptive error message
     */
    GeneralException(int errorCode, std::string message=""):errorCode(errorCode),message(message) {}

    virtual ~GeneralException(){}

    /**
     * Gets the error code
     * @return the error code
     */
    virtual int getErrorCode() const {
        return errorCode;
    }

    /**
     * Gets the message
     * @return the message
     */
    virtual const std::string& getMessage() const {
        return message;
    }
    
    /**
     * Sets the message
     * @param message the message
     */
    virtual void setMessage(const std::string& message) {
        this->message = message;
    }
    
    private:
    int errorCode;
    std::string message;
};

/* KB2 parser exceptions */

/**
 * Special exception for parsing errors.
 * The error code is always 1 (for parse error).
 * @author Bjoern Rabenstein
 */
class ParseException : public GeneralException {
public:
    /**
     * The constructor takes a line number and a file name as arguments to show where the parse error has occured.
     * The error code is always 1 (for parse error).
     * To the message saved in the exception the line number and file name are added.
     * @param lineNo the line number where the error has occured
     * @param fileName the name of the file where the error has occured
     * @param message the error message
     */
    ParseException(int lineNo, const std::string& fileName, const std::string& message):GeneralException(1,message) {
        std::ostringstream ost;
        ost << "parse error in line " << lineNo << " of " << fileName << ": " << GeneralException::getMessage();
        setMessage(ost.str());
    }

    virtual ~ParseException() {}
}
;

/**
 * Special exception if the parser is asked for a non-existent InputLine.
 * If this exception is thrown, the reason is usually a programming bug.
 * The error code is always 4.
 * @author Bjoern Rabenstein
 */
class NoSuchInputLineException : public GeneralException {
    public:
    /**
     * @param message explanatory error message - it is automatically appended by the string
     * <code> (This is usually a programming error.)</code>
     */
        NoSuchInputLineException(const std::string& message):GeneralException(4) {
            std::ostringstream ost;
            ost << message << " (This is usually a programming error.)";
            setMessage(ost.str());
        }

        virtual ~NoSuchInputLineException() {}
}
;

/**
 * Special exception for errors during construction of the prototype objects (Argument and InputLine) of the parser.
 * The error code is always 3 (for parser setup error).
 * @author Bjoern Rabenstein
 */
class ParserSetupException : public GeneralException {
public:
    /**
     * @param message the error message
     */
    ParserSetupException(const std::string& message):GeneralException(3) {
        setMessage(std::string("Error during setup of parser: ")+message);
    }

    virtual ~ParserSetupException() {}
};

/**
 * Special exception for type conversion errors.
 * The error code is always 2 (for type conversion error).
 * @author Bjoern Rabenstein
 */
class BadTypeConversionException : public GeneralException {
public:
    /**
     * @param originalType the type that was tried to be converted
     * @param targetType the requested type (to which the conversion was tried)
     */
    BadTypeConversionException(const std::string& originalType, const std::string& targetType):GeneralException(2) {
        std::ostringstream ost;
        ost << "you requested a variable of type " << targetType << ", but I cannot convert from " << originalType;
        setMessage(ost.str());
    }

    virtual ~BadTypeConversionException() {}
};

/* own exceptions */

/**
 * This exception is thrown if a file could not be opened.
 * @author Gernot Kieseritzky
 */

class FileIOException : public GeneralException 
{
public:
    /**
     * Constructor, error code is 6
     * @param filename Name file that had IO errors.
     * @param where Method that threw this Exception.
     */
    FileIOException(const string& filename, const string& where):GeneralException(10) 
	 {
	  setMessage("Could not open file "+filename+" in "+where+"!");
	 }

    virtual ~FileIOException(){}
};

/**
 * Simple Out of Range exception.
 * Thrown if tried to access past the end of a vector.
 * @author Gernot Kieseritzky
 */

class OutOfRangeException : public GeneralException 
{
public:
    /**
     * Constructor that takes a string of the method name that produced the error.
	  * Error code is 11
     * @param where Locus where exception was thrown.
     */
    OutOfRangeException(const string& where):GeneralException(11) 
	 {
	  setMessage("Out of Range Exception in "+where);
	 }

    virtual ~OutOfRangeException(){}
};

/**
 * Thrown if tried to access a non-existent state in instance of Titratable.
 * Error code is 12.
 * @author Gernot Kieseritzky
 */

class UnknownStateException : public GeneralException 
{
public:
    /**
     * Constructor.
     * @param state Non-existent state index.
	  * @param residue Titratable residue name.
	  * @param where String of method name where error occurred.
     */

    UnknownStateException(int state, const string& residue, const string& where):GeneralException(12)
	 {
 	  ostringstream ostr;
	  ostr << "Unknown state #" << state << (residue.empty()?"":" in ") << residue << " by " << where;

	  setMessage(ostr.str());
	 }

    virtual ~UnknownStateException(){}
};

/**
 * This exception is thrown if duplicate coordinates were found for an atom.
 * @author Gernot Kieseritzky
 */

class DuplicateCoordinatesException : public GeneralException 
{
public:
    /**
     * Constructor, error code is 13
     * @param atomno1 atom number 1
	  * @param atomno2 atom number 2
     * @param where Method that threw this Exception.
     */
    DuplicateCoordinatesException(int atomno1, int atomno2, const string& where):GeneralException(13)
	 {
	  ostringstream ostr;
	  ostr << "Duplicate coordinates for atoms with number " << atomno1 << " and " << atomno2 << " in " << where << ".";
	  setMessage(ostr.str());
	 }

    virtual ~DuplicateCoordinatesException(){}
};

/**
 * This exception is thrown if an atom or residue could not be found.
 * @author Gernot Kieseritzky
 * @ingroup Structure
 */

class AtomNotFoundException : public GeneralException 
{
public:
    /**
     * Constructor 1, error code is 14
     * @param what 'resid' or 'atomid'
	  * @param which illegal atom/residue number
     * @param where Method that threw this Exception.
     */
    AtomNotFoundException(const string& what, int which, const string& where):GeneralException(14)
	 {
	  ostringstream ostr;
	  ostr << "Unknown " << what << " with number " << which << " by " << where << ".";
	  setMessage(ostr.str());
	 }

    /**
    * Constructor 2
    * @param which illegal atom/residue name
    */
    AtomNotFoundException(const string& what, const string& which, const string& where):GeneralException(14)
	 {
	  ostringstream ostr;
	  ostr << "Unknown " << what << " named " << which << " by " << where << ".";
	  setMessage(ostr.str());
	 }

    virtual ~AtomNotFoundException(){}
};

/**
 * This exception is thrown if an unknown coordinate was tried to be evaluated.
 * @ingroup Potential
 * @ingroup Structure
 * @author Gernot Kieseritzky
 */

class PointNotFoundException : public GeneralException 
{
public:
    /**
     * Constructor, error code is 15
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param where Method that threw this exception.
     */
    PointNotFoundException(double x, double y, double z, const string& where) : GeneralException(15)
	 {
	  ostringstream ostr;
	  ostr << "Unknown point " << "(" << x << "," << y << "," << z << ") in " << where << ".";
	  setMessage(ostr.str());
	 }

    virtual ~PointNotFoundException(){}
};

/**
 * Exception that is thrown if there are ill-defined APBS parameters.
 * @author Gernot Kieseritzky
 */

class ParameterInconsistencyException : public GeneralException
{
public:
 /**
  * Constructor that takes a string with the error description.
  * @param explanation Error description.
  */
 ParameterInconsistencyException(const string& explanation):GeneralException(16) 
 {
  setMessage(explanation);
 }

 virtual ~ParameterInconsistencyException(){}
};

/**
 * Exception that is thrown if there is a problem in APBS.
 * @author Gernot Kieseritzky
 */

class APBSErrorException : public GeneralException
{
public:
 /**
  * Constructor that takes a description of the error and the APBS error code
  * @param apbsroutine APBS routine that failed
  * @param where method that threw exception.
  * @param APBSErrorCode APBS error code.
  */
 APBSErrorException(const string& apbsroutine, const string& where, int APBSErrorCode):GeneralException(17) 
 {
  ostringstream ostr;
  ostr << apbsroutine << " failed in " << where << ". APBS error code: " << APBSErrorCode;
  setMessage(ostr.str());
 }

 virtual ~APBSErrorException(){}
};

/**
 * Tried to run APBS without calling setup.
 * @author Gernot Kieseritzky
 */

class APBSWrapperSetupNotCalledException : public GeneralException
{
public:
 /**
  * Constructor that takes a description of the error and the APBS error code
  * @param where method that threw exception.
  */
 APBSWrapperSetupNotCalledException(const string& where):GeneralException(18) 
 {
  ostringstream ostr;
  ostr << where << ": Must call setup first!";
  setMessage(ostr.str());
 }

 virtual ~APBSWrapperSetupNotCalledException(){}
};

/**
 * Called setup twice.
 * @author Gernot Kieseritzky
 */

class APBSWrapperSetupCalledTwiceException : public GeneralException
{
public:
 /**
  * Constructor that takes a description of the error and the APBS error code
  */
 APBSWrapperSetupCalledTwiceException():GeneralException(19, "APBS_base::setup called twice!") {}

 virtual ~APBSWrapperSetupCalledTwiceException(){}
};

/**
 * Operator overloading to enable output of exceptions using operator "<<".
 * This implementation writes the result of GeneralException::getMessage and
 * of GeneralException::getErrorCode to the ostream.
 * @param o the output stream
 * @param e the exception
 * @return a reference to the output stream
 * author: Bjoern Rabenstein
 */
inline ostream& operator<<(ostream& o,GeneralException& e)
{
 o << e.getMessage() << " (Error " << e.getErrorCode() << ")";
 return o;
}
