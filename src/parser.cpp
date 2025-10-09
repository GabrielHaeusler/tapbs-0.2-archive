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
#include "parser.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdio.h>

std::vector<std::string> tokenize(const std::string& line){
    std::vector<std::string> result;
    std::istringstream ist(line);
    std::string token;
    while (ist>>token){
        result.push_back(token);
    }
    return result;
}

std::string getFirstToken(const std::string& line){
    std::string result;
    std::istringstream ist(line);
    if (ist>>result){
        if (result[0]!='#') return result;
    }
    return std::string();
}

IntegerArgument::IntegerArgument(bool required, long defaultValue, long lowerBound, long upperBound)
    : Argument(required),value(defaultValue),lowerBound(lowerBound),upperBound(upperBound) {
    if (!isRequired() && (value<lowerBound || value>upperBound)) {
        std::ostringstream ost;
        ost << "default value of integer argument is out of allowed range, default value is " << value <<
        ", upper Bound is " << upperBound << ", lower bound is " << lowerBound;
        throw ParserSetupException(ost.str());
    }
}

void IntegerArgument::parse(const std::string& text, int lineNo, const std::string& fileName) {
    long result;
    if(sscanf(text.c_str(),"%ld",&result)!=1){
        throw ParseException(lineNo, fileName,"argument is not an integer");
    }
    if (result < lowerBound){
        std::ostringstream ost;
        ost << "value of integer argument is to small, lower bound is " << lowerBound;
        throw ParseException(lineNo, fileName, ost.str());
    }
    if (result > upperBound){
        std::ostringstream ost;
        ost << "value of integer argument is to large, upper bound is " << upperBound;
        throw ParseException(lineNo, fileName, ost.str());
    }
    value = result;
}

FloatArgument::FloatArgument(bool required, double defaultValue, double lowerBound, double upperBound)
        : Argument(required),value(defaultValue),lowerBound(lowerBound),upperBound(upperBound) {
    if (!isRequired() && (value<lowerBound || value>upperBound)) {
        std::ostringstream ost;
        ost << "default value of float argument is out of allowed range, default value is " << value <<
        ", upper Bound is " << upperBound << ", lower bound is " << lowerBound;
        throw ParserSetupException(ost.str());
    }
}

void FloatArgument::parse(const std::string& text, int lineNo, const std::string& fileName) {
    double result;
    if(sscanf(text.c_str(),"%lf",&result)!=1){
        throw ParseException(lineNo, fileName,"argument is not a floating point number");
    }
    if (result < lowerBound){
        std::ostringstream ost;
        ost << "value of float argument is to small, lower bound is " << lowerBound;
        throw ParseException(lineNo, fileName, ost.str());
    }
    if (result > upperBound){
        std::ostringstream ost;
        ost << "value of float argument is to large, upper bound is " << upperBound;
        throw ParseException(lineNo, fileName, ost.str());
    }
    value = result;
}

void EnergyUnitArgument::parse(const std::string& text, int lineNo, const std::string& fileName) {
    std::string oldUnit = getStringValue();
    StringArgument::parse(text, lineNo, fileName);
    try {
        determineFactor();
    } catch (ParserSetupException e){
        StringArgument::parse(oldUnit, lineNo, fileName);
        throw ParseException(lineNo, fileName, e.getMessage());
    }
}

void EnergyUnitArgument::determineFactor() {
    std::string unit = getStringValue();
    if (unit=="" || unit=="kJ" || unit=="kJ/mol"){
        factor = 1.;
    } else if (unit=="kcal/mol" || unit=="kcal"){
        factor = KCAL;
    } else if (unit=="mV" || unit=="meV"){
        factor = MVOLT;
    } else if (unit=="V" || unit=="eV"){
        factor = VOLT;
    } else if (unit=="e^2/A" || unit=="e2A"){
        factor = E2A;
    } else if (unit=="pK"){
        factor = 1/LN10/KB;
    } else {
        throw ParserSetupException(std::string("unknown energy unit ")+unit);
    }
}

void InputLine::parse(const std::string& line, int lineNo, const std::string& fileName){
    if (line==""){
        if (defaultValue!="") parse(defaultValue, lineNo, fileName);
    } else {
        this->lineNo = lineNo;
        std::vector<std::string> tokens = tokenize(line);
        if (tokens.at(0)==keyword){
            if (tokens.size()-1 > args.size()){
                throw ParseException(lineNo, fileName, "too many arguments");
            }
            std::vector<std::string>::size_type i;
            for(i=1; i<tokens.size(); ++i){
                getArgument(i).parse(tokens[i],lineNo,fileName);
            }
            for ( ; i<=args.size();++i){
                if (getArgument(i).isRequired()){
                    throw ParseException(lineNo, fileName, "required argument missing");
                }
            }
        }
    }
}

InputParser::InputParser(const InputDescription& format, std::istream& in){
    using namespace std;
    currentUpperBound = input.end();
    currentLowerBound = input.end();
    cout << "I'm parsing stdin.\n";
    cout << "Your input is the following:\n\n";
    string line;
    int i = 0;
    while(getline(in,line)){
        cout << setw(4) << (++i) << ": " << line << std::endl;
        string firstToken = getFirstToken(line);
        if (firstToken != ""){
            // we have found something!
            InputDescription::const_iterator prototype = format.find(firstToken);
            if (prototype==format.end()){
                ostringstream ost;
                ost <<"unknown keyword '" << firstToken << "'";
                throw ParseException(i,"STDIN",ost.str());
            }
            if (input.find(firstToken)!=input.end() && !prototype->second.isUsableMoreThanOnce()){
                cout <<"WARNING: keyword '"<< firstToken << "' was already used in input, this line will be ignored.\n";
            }
            // now copy prototype:
            InputLine inputLine = prototype->second;
            inputLine.parse(line,i,"STDIN");
            input.insert(make_pair(firstToken,inputLine));
        }
    }
    // finally test for required keywords and add defaults for missing ones:
    for (InputDescription::const_iterator iter = format.begin(); iter != format.end(); ++iter){
        if (input.find(iter->first)==input.end()){
            // keyword is not yet in input
            const string& defaultLine = iter->second.getDefaultValue();
            if (defaultLine == "required"){
                ostringstream ost;
                ost << "required keyword '" << iter->first << "' not found in input.\n";
                throw ParseException(0,"STDIN",ost.str());
            }
            if (defaultLine != ""){
                InputLine inputLine = iter->second;
                inputLine.parse("",0,"STDIN");
                input.insert(make_pair(iter->first,inputLine));
            }
        }
    }
    cout << "\nI've finished parsing stdin.\n";
}

const InputLine& InputParser::getInputLine(const std::string& keyword) const {
    using namespace std;
    multimap<string,InputLine>::const_iterator iter = input.find(keyword);
    if (iter==input.end()){
        ostringstream ost;
        ost << "no input line available for keyword " << keyword;
        throw NoSuchInputLineException(ost.str());
    }
    return iter->second;
}

const InputLine& InputParser::getNextInputLine(){
    using namespace std;
    if (currentLowerBound==currentUpperBound){
        ostringstream ost;
        if (keyword.length()==0){
            ost << "no keyword set";
        } else {
            ost << "no more input lines available for keyword " << keyword;
        }
        throw NoSuchInputLineException(ost.str());
    }
    return (currentLowerBound++)->second;
}

double InputLine::getEnergy(int firstArgument, double temp) const{
    double interaction = getArgument(firstArgument).getDoubleValue();
    return interaction/(dynamic_cast<EnergyUnitArgument&>(getArgument(firstArgument+1))).getFactor(temp);
}
