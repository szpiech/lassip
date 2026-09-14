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
#include "param_t.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

bool param_t::addFlag(string flag, bool value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        if (value) buffer = "true";
        else buffer = "false";
        argb[flag] = value;
        help[flag] = "<bool>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, double value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%.2f", value);
        buffer = charBuffer;
        argd[flag] = value;
        help[flag] = "<double>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%d", value);
        buffer = charBuffer;
        argi[flag] = value;
        help[flag] = "<int>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, char value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%c", value);
        buffer = charBuffer;
        argch[flag] = value;
        help[flag] = "<char>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        args[flag] = value;
        help[flag] = "<string>: " + description + "\n\tDefault: " + value;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, const char value[], string label, string description)
{
    return this->addFlag(flag, string(value), label, description);
}

bool param_t::addListFlag(string flag, double value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%f", value);
        buffer = charBuffer;
        listargd[flag].push_back(value);
        help[flag] = "<double1> ... <doubleN>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}
bool param_t::addListFlag(string flag, char value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%c", value);
        buffer = charBuffer;
        listargch[flag].push_back(value);
        help[flag] = "<char1> ... <charN>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addListFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        snprintf(charBuffer, sizeof(charBuffer), "%d", value);
        buffer = charBuffer;
        listargi[flag].push_back(value);
        help[flag] = "<int1> ... <intN>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool  param_t::addListFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        listargs[flag].push_back(value);
        help[flag] = "<string1> ... <stringN>: " + description + "\n\tDefault: " + value;
        labels[flag] = label;
        flagOrder.push_back(flag);
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addListFlag(string flag, const char value[], string label, string description)
{
    return this->addListFlag(flag, string(value), label, description);
}

void param_t::printHelp()
{
    //Help goes to stdout: it is what the user asked for, not a diagnostic, and
    //this lets `lassip --help | less` work.
    if (usage.length() > 0) cout << usage << "\n";
    cout << preamble << endl;

    //Flags are grouped by the label given to addFlag, in the order the groups
    //were first registered, so that related options appear together instead of
    //alphabetically -- which used to put --avg-spec first and scatter the two
    //stages of the workflow through one long list.
    vector<string> groupOrder;
    map<string, vector<string> > byGroup;
    for (unsigned int i = 0; i < flagOrder.size(); i++)
    {
        string flag = flagOrder[i];
        if (help.count(flag) == 0) continue;
        string group = labels[flag];
        if (group.compare("SILENT") == 0) continue;
        if (byGroup.count(group) == 0) groupOrder.push_back(group);
        byGroup[group].push_back(flag);
    }

    for (unsigned int g = 0; g < groupOrder.size(); g++)
    {
        string group = groupOrder[g];
        if (group.length() > 0) cout << "----------" << group << "----------\n\n";
        else cout << "----------Command Line Arguments----------\n\n";
        for (unsigned int i = 0; i < byGroup[group].size(); i++)
        {
            string flag = byGroup[group][i];
            cout << flag << " " << help[flag] << "\n\n";
        }
    }

    return;
}

void param_t::setVersion(string str)
{
    version = str;
    return;
}

void param_t::setUsage(string str)
{
    usage = str;
    return;
}

bool param_t::isFlagSet(string flag)
{
    return (isSet.count(flag) > 0 && isSet[flag]);
}

bool param_t::goodDouble(string str)
{
    string::iterator it;
    //int dashCount = 0;
    int decimalCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '.' && *it != '-') return 0;
        if (*it == '.') decimalCount++;
        if (*it == '-' && it != str.begin()) return 0;
        if (/*dashCount > 1 || */decimalCount > 1) return 0;
    }
    return 1;
}

bool param_t::goodInt(string str)
{
    string::iterator it;
    //int dashCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '-') return 0;
        if (*it == '-' && it != str.begin()) return 0;
        //if (dashCount > 1) return 0;
    }
    return 1;
}

bool param_t::goodChar(string str)
{
    if (str.length() > 1) return 0;
    return 1;
}

bool param_t::parseCommandLine(int argc, char *argv[])
{
    int badFlags = 0;

    for (int i = 1; i < argc; i++)
    {
        if (isSet.count(argv[i]) > 0)
        {
            cerr << "ERROR: Duplicate " << argv[i] << " found.\n";
            badFlags++;
            break;
        }
        else if (argb.count(argv[i]) > 0)
        {
            argb[argv[i]] = !argb[argv[i]];
            isSet[argv[i]] = true;
        }
        else if (argi.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodInt(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                badFlags++;
                break;
            }
            else
            {
                argi[argv[i]] = atoi(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargi.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargi[argv[i]].clear();//clear the default value
                isSet[argv[i]] = true;
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodInt(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargi[argv[flagIndex]].push_back(atoi(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodInt(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargi[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argd.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodDouble(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid double.\n";
                badFlags++;
                break;
            }
            else
            {
                argd[argv[i]] = atof(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargd.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargd[argv[i]].clear();//clear the default value
                isSet[argv[i]] = true;
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodDouble(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargd[argv[flagIndex]].push_back(atof(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodDouble(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid double.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargd[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (args.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                args[argv[i]] = argv[i + 1];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargs.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargs[argv[i]].clear();//clear the default value
                isSet[argv[i]] = true;
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (argv[i + 1][0] != '-') //make sure the next value isn't another flag
                    {
                        listargs[argv[flagIndex]].push_back(string(argv[i + 1]));
                        i++;
                    }
                    else //if the next value is another flag...
                    {
                        if (listargs[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argch.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodChar(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid character.\n";
                badFlags++;
                break;
            }
            else
            {
                argch[argv[i]] = argv[i + 1][0];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargch.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargch[argv[i]].clear();//clear the default value
                isSet[argv[i]] = true;
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodChar(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargch[argv[flagIndex]].push_back(atoi(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodChar(string(argv[i + 1])) && !flagExists(string(argv[i + 1])))
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid character.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargch[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (ARG_HELP_SHORT.compare(argv[i]) == 0)
        {
            argb[ARG_HELP] = true;
            isSet[ARG_HELP] = true;
        }
        else //if (argv[i][0] == '-')
        {
            cerr << "ERROR: command line flag " << argv[i] << " not recognized.\n";
            cerr << "Run lassip --help for a list of options.\n";
            badFlags++;
        }
    }

    if (getBoolFlag(ARG_HELP))
    {
        this->printHelp();
        throw ParamExit(0);
    }

    if (getBoolFlag(ARG_VERSION))
    {
        cout << version << "\n";
        throw ParamExit(0);
    }

    if (badFlags) throw ParamExit(EXIT_USAGE);

    return 0;
}

bool param_t::flagExists(string flag)
{
    return (help.count(flag) > 0);
}

param_t::param_t()
{
    this->addFlag(ARG_HELP, false, "SILENT", "Prints this help dialog.");
    this->addFlag(ARG_VERSION, false, "SILENT", "Prints the version and exits.");
}

bool param_t::getBoolFlag(string flag)
{
    if (argb.count(flag) > 0) return argb[flag];

    cerr << "ERROR: There are no bool flags named " << flag << "\n";
    throw 0;
}

double param_t::getDoubleFlag(string flag)
{
    if (argd.count(flag) > 0) return argd[flag];

    cerr << "ERROR: There are no double flags named " << flag << "\n";
    throw 0;
}

int param_t::getIntFlag(string flag)
{
    if (argi.count(flag) > 0) return argi[flag];

    cerr << "ERROR: There are no int flags named " << flag << "\n";
    throw 0;
}

char param_t::getCharFlag(string flag)
{
    if (argch.count(flag) > 0) return argch[flag];

    cerr << "ERROR: There are no char flags named " << flag << "\n";
    throw 0;
}

string param_t::getStringFlag(string flag)
{
    if (args.count(flag) > 0) return args[flag];

    cerr << "ERROR: There are no string flags named " << flag << "\n";
    throw 0;
}

vector<string> param_t::getStringListFlag(string flag)
{
    if (listargs.count(flag) > 0) return listargs[flag];

    cerr << "ERROR: There are no string list flags named " << flag << "\n";
    throw 0;
}

vector<int> param_t::getIntListFlag(string flag)
{
    if (listargi.count(flag) > 0) return listargi[flag];

    cerr << "ERROR: There are no int list flags named " << flag << "\n";
    throw 0;
}

vector<double> param_t::getDoubleListFlag(string flag)
{
    if (listargd.count(flag) > 0) return listargd[flag];

    cerr << "ERROR: There are no double list flags named " << flag << "\n";
    throw 0;
}

vector<char> param_t::getCharListFlag(string flag)
{
    if (listargch.count(flag) > 0) return listargch[flag];

    cerr << "ERROR: There are no int list flags named " << flag << "\n";
    throw 0;
}

void param_t::setPreamble(string str)
{
    preamble = str;
    return;
}
