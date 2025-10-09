/***************************************************************************
 *   Copyright (C) 2004 by Bjoern Rabenstein                               *
 *   rabe@chemie.fu-berlin.de                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef KARLSBERGPARSER_H
#define KARLSBERGPARSER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <istream>
#include <sstream>
#include <map>
#include <vector>
#include <limits>
#include <memory>

#include "utils.h"
#include "exceptions.h"

using namespace std;

/// k_B in kJ/(mol*K), actually R. From 1986 CRC Handbook.
const double KB = 0.008314510;

/// 1 kJ/mol in kcal/mol
const double KCAL = 0.239006;

/// 1 kJ/mol in eV/e or V
const double VOLT   = 0.01036427;

/// 1 kJ/mol in meV/e or mV
const double MVOLT  = 10.36427;

/// 1 kJ/mol in e^2/(A*mol)
const double E2A  = 7.197585e-4;

/**
 * @file parser.h contains various classes for parsing input files.
 */

class InputLine;

/**
 * This map describes the format of the input to stdin.
 * The keys of this map are the possible keywords, the values
 * are InputLine prototype objects describing the keyword and its possible arguments.
 */
typedef std::map<std::string, InputLine> InputDescription;

/**
 * Splits a line at whitespaces and puts all token into a std::vector.
 * @param line the line to split
 * @return the token vector, elemtns are of type std::string
 */
std::vector<std::string> tokenize(const std::string& line);

/**
 * Returns the first token (white-space seperated word) of the given string.
 * If the string is empty or the first token starts with <code>#</code>,
 * an empty string is returned.
 * @param line the line to parse
 * @return the first token or an empty string
 */
std::string getFirstToken(const std::string& line);

/**
 * Abstract base class describing an argument of a keyword. Concrete instances
 * of this class are used to represent arguments of the input and also
 * (as prototype objects) to describe the format of an InputLine.
 */
class Argument {
    bool required;
protected:
    /**
     * Protected constructor to be used by derived classes.
     * @param required is this argument required?
     */
    Argument(bool required=true)
            : required(required) {}
public:
    /**
     * Virtual copy constructor
     * @return a pointer to the clone
     */
    virtual Argument* clone() const=0;

    /**
     * Is this argument required?
     */
    virtual bool isRequired() const {
        return required;
    }

    virtual ~Argument() {}

    /**
     * Parses the given argument text, throwing an ParseException if it does not
     * follow the given rules. In case of success, the internal value of this argument
     * Object is set appropriately.
     * @param text the argument text to parse
     * @param lineNo the line number from where the text is taken (for generation of error messages)
     * @param fileName the file name from wehre the text is taken (for generation of error messages)
     */
    virtual void parse(const std::string& text, int lineNo, const std::string& fileName)=0;

    /**
     * Gets the value of this argument as a std::string.
     * @return the value of the argument as a string
     */
    virtual const std::string& getStringValue() const=0;

    /**
     * Gets the value of this argument as an <code>long</code>.
     * Might involve rounding. An exception is thrown for
     * arguments whose value cannot be expressed as <code>long</code>.
     * @return the value of the argument as <code>long</code>
     */
    virtual long getLongValue() const=0;

    /**
     * Gets the value of this argument as an <code>double</code>.
     * An exception is thrown for
     * arguments whose value cannot be expressed as <code>double</code>.
     * @return the value of the argument as <code>double</code>
     */
    virtual double getDoubleValue() const=0;
};

/**
 * An argument that is an arbitrary string, represented internally as std::string object.
 */
class StringArgument : public Argument {
    std::string value;
public:
    /**
    * This constructor is used to construct the protoype objects describing an argument.
    * @param required is this argument required?
     * @param defaultValue the default value of this argument, empty string means no default value
     */
    StringArgument(bool required=true, const std::string& defaultValue="")
            : Argument(required), value(defaultValue) {}

    virtual void parse(const std::string& text, int lineNo, const std::string& fileName) {
        value = text;
    }

    virtual const std::string& getStringValue() const {
        return value;
    }

    virtual Argument* clone() const {
        return new StringArgument(*this);
    }

    /**
     * Always throws an BadTypeConversionException.
     */
    virtual long getLongValue() const {
        throw BadTypeConversionException("std::string","int");
        return 0;
    }

    /**
     * Always throws an BadTypeConversionException.
     */
    virtual double getDoubleValue() const {
        throw BadTypeConversionException("std::string","double");
        return 0.;
    }
};

/**
 * An argument that is an integer number, represented internally as <code>long</code>.
 */
class IntegerArgument : public Argument {
    long value, lowerBound, upperBound;
    std::string stringValue;
public:
    /**
    * This constructor is used to construct the protoype objects describing an argument.
    * @param required is this argument required?
     * @param defaultValue the default value of this argument - if <code>numeric_limits&lt;int&gt;::min()</code>,
     * there is no default value. If <i>required</i> is <code>true</code>, this parameter is ignored.
     * @param lowerBound the minimum allowed value
     * @param upperBound the maximum allowed value.
     */
    IntegerArgument(bool required=true,
                    long defaultValue = std::numeric_limits<long>::min(),
                    long lowerBound = std::numeric_limits<long>::min(),
                    long upperBound = std::numeric_limits<long>::max());

    virtual Argument* clone() const {
        return new IntegerArgument(*this);
    }

    virtual void parse(const std::string& text, int lineNo, const std::string& fileName);

    virtual const std::string& getStringValue() const {
        return stringValue;
    }

    virtual long getLongValue() const {
        return value;
    }

    virtual double getDoubleValue() const {
        return value;
    }
};

/**
* An argument that is a floating point number, represented internally as <code>double</code>
*/
class FloatArgument : public Argument {
    double value, lowerBound, upperBound;
    std::string stringValue;
public:
    /**
    * This constructor is used to construct the protoype objects describing an argument.
    * @param required is this argument required?
    * @param defaultValue the default value of this argument - if <code>numeric_limits&lt;double&gt;::min()</code>,
    * there is no default value. If <i>required</i> is <code>true</code>, this parameter is ignored.
    * @param lowerBound the minimum allowed value
    * @param upperBound the maximum allowed value.
     */
    FloatArgument(bool required=true,
                  double defaultValue = -std::numeric_limits<double>::max(),
                  double lowerBound = -std::numeric_limits<double>::max(),
                  double upperBound = std::numeric_limits<double>::max());

    virtual Argument* clone() const {
        return new FloatArgument(*this);
    }

    virtual void parse(const std::string& text, int lineNo, const std::string& fileName);

    virtual const std::string& getStringValue() const {
        return stringValue;
    }

    virtual long getLongValue() const {
        return static_cast<long>(value);
    }

    virtual double getDoubleValue() const {
        return value;
    }

};

/**
 * A special StringArgument, the value of the argument must be one of these:
 * <ul>
    * <li><code>kJ/mol</code> or <code>kJ</code> for kJ/mol </li>
    * <li><code>kcal/mol</code> or  <code>kcal</code> for kcal/mol</li>
    * <li><code>mV</code> or <code>meV</code> for mV or meV</li>
    * <li><code>V</code> or <code>eV</code> for V or eV</li>
    * <li><code>e^2/A</code> or <code>e2A</code> for e<sup>2</sup>/&Aring;</li>
 * <li><code>pK</code> for pK units</li>
    * </ul>
 */
class EnergyUnitArgument : public StringArgument {
    double factor;
public:
    /**
     * This constructor is used to construct the protoype objects describing an argument.
     * @param required is this argument required?
     * @param defaultValue the default value of this argument, empty string means no default value - must
     * be an allowed value for this type of argument.
     */
    EnergyUnitArgument(bool required=true, const std::string& defaultValue="")
            : StringArgument(required, defaultValue) {
        determineFactor();
    }
    virtual Argument* clone() const {
        return new EnergyUnitArgument(*this);
    }

    virtual void parse(const std::string& text, int lineNo, const std::string& fileName);
    /**
     * Gets the factor for the energy unit represented by this argument. The factor is
     * the amount of the energy 1 kJ/mol in the given energy unit.
     * @param temperature For some energy units, the conversion to/from kJ/mol needs
     * a temperature in K, which is given by the value of this parameter. For unit conversions
     * where the temperature is not used, this parameter has no effect.
     * @return the factor for this energy unit.
     */
    virtual double getFactor(double temperature) const {
        return (getStringValue()=="pK"? factor/temperature : factor);
    }

private:
    void determineFactor();
};

/**
 * This class represents an input line, i.e. a keyword and its arguments.
 * It is also used as a prototype object to describe the format of an input line.
 */
class InputLine {
    std::string keyword;
    std::vector<Argument*> args;
    std::string defaultValue;
    bool usableMoreThanOnce;
    int lineNo;
public:
    /**
     * This constructor is used to construct the protoype objects describing the input line format.
     * @param keyword the keyword
     * @param defaultValue If the keyword of this InputLine is not present in the input, the string
     * given as this defaultValue is taken instead as a complete input line to be parsed by this InputLine object.
     * However, there are a few special cases:<ul>
     * <li>If the defaultValue is an empty string (and the keyword of this InputLine is not present in the input),
     * no parsing will be done. A keyword not in the input
     * will simply be not there.</li>
     * <li>If the defaultValue is <code>required</code>
     * (and the keyword of this InputLine is not present in the input),
     * no parsing will be done either. Instead, an exception will be thrown.</li>
     * @param usableMoreThanOnce if true, the described keyword may appear more than once in the input. If false,
     * only the first appearance will be recognized. Other appearances will be ignored (but a warning will be issued).
     */
    InputLine(const std::string& keyword,
              const std::string& defaultValue="",
              bool usableMoreThanOnce=false)
            :keyword(keyword), defaultValue(defaultValue), usableMoreThanOnce(usableMoreThanOnce) {}

    /**
     * Default constructor
    */
    InputLine() {}

    /**
     * Copy constructor
     */
    InputLine(const InputLine& that) {
        copyMembers(that);
    }

    /**
     * Destructor deletes the arguments.
     */
    ~InputLine() {
        for(std::vector<Argument*>::iterator iter = args.begin(); iter != args.end(); ++iter) {
            delete(*iter);
        }
    }

    /**
     * Assignment operator
     */
    InputLine& operator=(const InputLine& that) {
        if (this != &that) {
            for(std::vector<Argument*>::iterator iter = args.begin(); iter != args.end(); ++iter) {
                delete(*iter);
            }
            args.clear();
            copyMembers(that);
        }
        return *this;
    }

    /** Returns the number of the input line from where this InputLine was created. */
    int getLineNo() const {
        return lineNo;
    }
    
    /**
     * Adds an argument to this InputLine. It is used to assemble the prototype InputLine object.
     * The Argument is <em>cloned</em> by this method.
     * @param arg the Argument
     * @return a reference to this InputLine.
     */
    InputLine& addArgument(const Argument& arg) {
        args.push_back(arg.clone());
        return *this;
    }

    /**
     * Gets the Argument with the given (1-based) index.
     * @param index the 1-based index of the argument to return
     * @return a reference to the Argument
     */
    Argument& getArgument(int index) const {
        return *(args.at(index-1));
    }

    /**
     * Gets the keyword of this input line.
     * @return the keyword
     */
    const std::string& getKeyword() const {
        return keyword;
    }

    /**
     * Gets the default value for the whole input line.
     */
    const std::string& getDefaultValue() const {
        return defaultValue;
    }
    
    /**
     * Is the keyword of this InputLine usable more than once?
     */
    bool isUsableMoreThanOnce() const {
        return usableMoreThanOnce;
    }

    /**
     * Adds this InputLine to the given InputDescription map.
     * @param inputDescription the InputDescription where to add ourselves
     */
    void addToInputDescription(InputDescription& inputDescription) const {
        inputDescription[keyword] = *this;
    }

    /**
     * Parses an inpur line. The values of the Argument objects are set
     * according to the parse result.
     * @param line The line of input to parse. If the line is an empty string, the default
     * value for the input line is parsed instead. If the first item of the line is not equal to
     * the keyword, no parsing will occur.
     * @param lineNo the line number from where the text is taken (for generation of error messages)
     * @param fileName the file name from wehre the text is taken (for generation of error messages)
     */
    void parse(const std::string& line, int lineNo, const std::string& fileName);

    /**
     * Convenieance method to extract an energy value (in kJ/mol)
     * from two consecutive arguments, first the double value of the energy,
     * second the unit.
     * @param firstArgument the first argument (the energy amount as double) - the following argument is supposed to
     * be the energy unit
     * @param temp the temperature used to convert pK units
     * @return the energy value in kJ/mol
     */
    double getEnergy(int firstArgument, double temp) const;

private:
    void copyMembers(const InputLine& that) {
        this->defaultValue = that.defaultValue;
        this->keyword = that.keyword;
        this->usableMoreThanOnce = that.usableMoreThanOnce;
        this->lineNo = that.lineNo;
        for(std::vector<Argument*>::const_iterator iter = that.args.begin(); iter != that.args.end(); ++iter) {
            this->addArgument(**iter);
        }
    }

};

/**
 * The parser for the input to stdin.
 */
class InputParser {
    /** keys of this map are keywords, values are InputLine objects */
    std::multimap<std::string,InputLine> input;
    std::multimap<std::string,InputLine>::const_iterator currentLowerBound;
    std::multimap<std::string,InputLine>::const_iterator currentUpperBound;
    std::string keyword;
public:
    /**
     * The constructor takes the input stream to read from (usually std::cin) and an karlsberg::InputFormat that
     * describes the format of the input. The input stream is read and all parsing is done during construction
     * of the parser
     * @param format the input format description
     * @param in the input stream to read from
     */
    InputParser(const InputDescription& format, std::istream& in);

    /**
     * Returns true if there is at least one karlsberg::InputLine present for the given keyword.
     * This might also be the default line.
     */
    bool hasInputLine(const std::string& keyword) const {
        return input.count(keyword)>0;
    }
    
    /**
     * Gets the first karlsberg::InputLine for the given keyword.
     * @throw karlsberg::NoSuchInputLineException if no karlsberg::InputLine
     * can be found for the given keyword
     */
    const InputLine& getInputLine(const std::string& keyword) const;

    /**
     * Sets a keyword to be relevant for InputParser::hasNextInputLine and InputParser::getNextInputLine.
     * Setting the keyword resets the iterator used by these two methods.
     */
    void setKeyword(const std::string& keyword){
        using namespace std;
        this->keyword = keyword;
        pair<multimap<string,InputLine>::const_iterator,multimap<string,InputLine>::const_iterator>
                range = input.equal_range(keyword);
        currentLowerBound = range.first;
        currentUpperBound = range.second;
    }

    /**
     * Returns true if a keyword has been set before using InputParser::setKeyword and the next call of
     * InputParser::getNextInputLine will succeed, if InputParser::setKeyword is not called again.
     */
    bool hasNextInputLine() const {
        return currentLowerBound != currentUpperBound;
    }

    /**
     * Gets an InputLine suitable for the keyword set before using InputParser::setKeyword.
     * The first call of this method after calling InputParser::setKeyword returns the first suitable
     * InputLine. Succesive calls return the following InputLines for the same keyword, if any.
     * Thus, this method together with InputParser::setKeyword and InputParser::hasNextInputLine
     * can be used to retrieve multiple occurrences of InputLines for a specific keyword.
     * @throw karlsberg::NoSuchInputLineException if no keyword was set or no further
     * InputLine is available
     */
    const InputLine& getNextInputLine();
    
};

#endif
