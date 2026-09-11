/* param_t -- a class for basic command line argument parsing
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#ifndef __PARAM_T_H__
#define __PARAM_T_H__

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cctype>

using namespace std;

//sysexits.h conventions, so that a caller can tell a bad command line from bad
//data from an I/O failure.
const int EXIT_USAGE   = 64;
const int EXIT_DATAERR = 65;
const int EXIT_IOERR   = 74;
const int EXIT_INTERNAL = 70;

const string ARG_HELP = "--help";
const string ARG_HELP_SHORT = "-h";
const string ARG_VERSION = "--version";

//Thrown by parseCommandLine when the program should stop before doing any work:
//code 0 after --help or --version, EXIT_USAGE after a command line error.
struct ParamExit
{
    int code;
    ParamExit(int c) : code(c) {}
};

class param_t
{
public:

    /*
      The addFlag and addListFlag methods initialize a 
    */
    bool addFlag(string flag, bool value, string label, string description);
    bool addFlag(string flag, double value, string label, string description);
    bool addFlag(string flag, int value, string label, string description);
    bool addFlag(string flag, char value, string label, string description);
    bool addFlag(string flag, string value, string label, string description);
    bool addFlag(string flag, const char value[], string label, string description);

    bool addListFlag(string flag, string value, string label, string description);
    bool addListFlag(string flag, const char value[], string label, string description);
    bool addListFlag(string flag, int value, string label, string description);
    bool addListFlag(string flag, double value, string label, string description);
    bool addListFlag(string flag, char value, string label, string description);

    void printHelp();
    void setVersion(string str);
    void setUsage(string str);

    //true if the flag was given on the command line, as opposed to holding its
    //default. Lets callers distinguish "not supplied" from "supplied the
    //default value" without comparing against sentinel strings.
    bool isFlagSet(string flag);

    bool parseCommandLine(int argc, char *argv[]);

    bool getBoolFlag(string flag);
    double getDoubleFlag(string flag);
    int getIntFlag(string flag);
    char getCharFlag(string flag);
    string getStringFlag(string flag);

    vector<string> getStringListFlag(string flag);
    vector<int> getIntListFlag(string flag);
    vector<double> getDoubleListFlag(string flag);
    vector<char> getCharListFlag(string flag);

    void setPreamble(string str);

    param_t();


private:

    map<string, bool> argb;
    map<string, double> argd;
    map<string, int> argi;
    map<string, char> argch;
    map<string, string> args;

    map<string, vector< string > > listargs;
    map<string, vector< int > > listargi;
    map<string, vector< double > > listargd;
    map<string, vector< char > > listargch;

    map<string, string> help;
    map<string, bool> isSet;
    map<string, string> labels;

    bool goodDouble(string str);
    bool goodInt(string str);
    bool goodChar(string str);

    bool flagExists(string flag);

    string preamble;
    string version;
    string usage;
    vector<string> flagOrder;
};

#endif
