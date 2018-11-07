// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#include <algorithm> 
#include <cctype>
#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>

#include "MPar.h"

#include "MParManager.h"
#include "MParContainer.h"
#include "MPar.h"

/** \class MParManager
\ingroup lib_base
Parameters Manager class. Stores and dumps all parameters used by the framework.
It's a singleton class of which object can be obtained using instance() method.
Paramaters mamager must be initializatied before
MDetectorManager::initParameterContainers() is called since it checks whether
the requested parameter containers exists.
*/

// for trim functions see
// https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring

/** Trim from start (in place)
 *
 * \param s string
 */
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

/** Trim from end (in place)
 *
 * \param s string
 */
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

/** Trim from both ends (in place)
 *
 * \param s string
 */
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

/** Trim from start (copying)
 *
 * \param s string
 * \return trimmed string
 */
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

/** Trim from end (copying)
 *
 * \param s string
 * \return trimmed string
 */
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

/** Trim from both ends (copying)
 *
 * \param s string
 * \return trimmed string
 */
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

/** Remove all tabs (in place)
 *
 * \param s string
 */
void simplify(std::string & s)
{
    size_t pos = 0;
    while(1)
    {
        pos = s.find_first_of('\t', pos);
        if (pos == s.npos)
            return;
        s.replace(pos, 1, " ");
    }
}

/** Check if float.
 *
 * \sa Stackoverflow #447206
 *
 * \param str string
 * \return is float
 */
bool isFloat(const string & str)
{
    std::istringstream iss(str);
    float f;
    iss >> noskipws >> f;
    return iss.eof() && !iss.fail();
}

MParManager * MParManager::pm = nullptr;

/** Returns instance of the Detector Manager class.
 *
 * \return manager instance
 */
MParManager * MParManager::instance()
{
    if (!pm)
        pm = new MParManager;

    return pm;
}

/** Shortcut
 * \return MParManager instance
 */
MParManager * pm()
{
    return MParManager::instance();
}

/** Default constructor
 */
MParManager::MParManager()
{
}

/** Destructor
 */
MParManager::~MParManager()
{
}

/** Parse source file
 *
 * \return success
 */
bool MParManager::parseSource()
{
    ifstream ifs(source);
    size_t length = 0;
    if (ifs)
    {
        ifs.seekg(0, ifs.end);
        length = ifs.tellg();
        ifs.seekg(0, ifs.beg);
    }

    char * cbuff = new char[length];

    WhatNext wn = WNContainer;

    std::string cont_name;
    std::string param_name;
    std::string type_name;
    std::vector<std::string> values;

    while(!ifs.eof())
    {
        ifs.getline(cbuff, length);

        std::string str(cbuff);
        trim(str);
        simplify(str);

        size_t pos = 0, pos2 = 0;

        // check if comment or empty line
        if (str[0] == '#' or (str.length() == 0 and wn != WNParamCont))
        {
            continue;
        }
        // if container mark found, check whether it should be there, e.g. container after param new line is forbidden
        else if (str[0] == '[')
        {
            if (wn == WNContainer or wn == WNContainerOrParam)
            {
                pos = str.find_first_of(']', 1);
                cont_name = str.substr(1, pos-1);
//                 printf("Found container %s\n", cont_name.c_str());

                containers.insert(std::pair<std::string, MParContainer *>(cont_name, new MParContainer(cont_name)));

                param_name.clear();

                wn = WNContainerOrParam;
                continue;
            }
            else
            {
                std::cerr << "Didn't expected container here: " << std::endl << str << std::endl;
                return false;
            }
        }
        else
        {
            // check if container name
            if (wn == WNContainer)
            {
                std::cerr << "Expected container name here: " << std::endl << str << std::endl;
                return false;
            }
            else if (wn == WNParam or wn == WNContainerOrParam)
            {
                // find parameter name ended with :
                pos = str.find_first_of(':', 1);
                param_name = str.substr(0, pos);
                trim(param_name);

                // prepare other variables
                type_name.clear();
                values.clear();

                // find type name
                pos = str.find_first_not_of(' ', pos+1);

                // type name must be in the same line like param name
                if (str[pos] == '\\')
                {
                    std::cerr << "No type name detected in the param name line:" << std::endl << str << std::endl;
                    return false;
                }
                pos2 = str.find_first_of(' ', pos+1);

                type_name = str.substr(pos, pos2-pos);
                if (type_name != "Int_t" and type_name != "Float_t" and type_name != "Double_t")
                {
                    std::cerr << "Invalid param type '" << type_name << "' in line:" << std::endl << str << std::endl;
                    return false;
                }

//                 printf("Found param=%s with type=%s\n", param_name.c_str(), type_name.c_str());

                wn = WNParamCont;
            }

            if (wn == WNParamCont)
            {
//                 printf("Search for values in %s\n", str.substr(pos2, -1).c_str());
                wn = parseValues(str.substr(pos2, -1), values);
                if (wn == WNContainerOrParam)
                {
                    if (values.size() == 0)
                    {
                        std::cerr << "No values given in line:" << std::endl << str << std::endl;
                        return false;
                    }
                    else
                    {
                        containers[cont_name]->initParam(param_name, type_name, values);
                    }
                }
            }
        }
    }

    return true;
}

/** Parse single value.
 *
 * \param str string with values
 * \param values vector to store values
 * \return next parsing step
 */
MParManager::WhatNext MParManager::parseValues(const std::string & str, std::vector<std::string> & values)
{
    size_t pos = 0, pos2 = 0;

    while(1)
    {
        pos = str.find_first_not_of(' ', pos2);
//         printf("Found %c at %d\n", str[pos], pos);
        // if new line detected
        if (str[pos] == '\\')
        {
            return WNParamCont;
        }

        // if end of line
        if (pos == str.npos)
        {
            break;
        }

        pos2 = str.find_first_of(' ', pos+1);
        std::string match = str.substr(pos, pos2-pos);
        if (isFloat(match))
        {
            values.push_back(match);
        }
        else
        {
            std::cerr << "Value is not a number:" << std::endl << str << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    return WNContainerOrParam;
}

/** Write params to destination.
 *
 * Needs to be implemented.
 */
void MParManager::writeDestination() const
{
}

/** Print containers
 */
void MParManager::print() const
{
     std::map<std::string, MParContainer *>::const_iterator it = containers.begin();
     for (; it != containers.end(); ++it)
         it->second->print();
}

/** Add new parameter container.
 *
 * \param cont_name container name
 * \param parcont container object
 * \return success
 */
bool MParManager::addParameterContainer(const std::string& cont_name, MPar* parcont)
{
    MParContainer * pc = containers[cont_name];
    if (!pc)
    {
        std::cerr << "Container " << cont_name << " doesn't exists!" << std::endl;
        exit(EXIT_FAILURE);
    }

    parcont->clear();
    bool ret = parcont->getParams(pc);
    if (!ret)
    {
        std::cerr << "Initialization of " << cont_name << " param container failed!" << std::endl;
        exit(EXIT_FAILURE);
    }

    parconts.insert(std::pair<std::string, MPar*>(cont_name, parcont));
    return true;
}

/** Get parameter container by name.
 *
 * \param cont_name container name
 * \return pointer to container
 */
MPar * MParManager::getParameterContainer(const std::string& cont_name)
{
    return parconts[cont_name];
}
