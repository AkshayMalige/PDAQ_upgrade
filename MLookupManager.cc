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


#include "MLookupManager.h"
#include "MLookupContainer.h"
#include "MLookup.h"

/** class MLookupManager
\ ingroup lib_base

Parameters Manager class. Stores and dumps all parameters used by the framework.

It's a singleton class of which object can be obtained using instance() method.

Paramaters mamager must be initializatied before
MDetectorManager::initParameterContainers() is called since it checks whether
the requested parameter containers exists.
*/
using namespace std;
/*
extern void trim(std::string &s);
extern void simplify(std::string & s);
extern bool isFloat(const string & str);*/


 void Ltrim(std::string &s);
 void Lsimplify(std::string & s);
 bool LlisFloat(const string & str);


/** Trim from start (in place)
 *
 * \param s string
 */
static inline void lLtrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

/** Trim from end (in place)
 *
 * \param s string
 */
static inline void rLtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

/** Trim from both ends (in place)
 *
 * \param s string
 */
inline void Ltrim(std::string &s) {
    lLtrim(s);
    rLtrim(s);
}

/** Trim from start (copying)
 *
 * \param s string
 * \return Ltrimmed string
 */
static inline std::string lLtrim_copy(std::string s) {
    lLtrim(s);
    return s;
}

/** Trim from end (copying)
 *
 * \param s string
 * \return trimmed string
 */
static inline std::string rLtrim_copy(std::string s) {
    rLtrim(s);
    return s;
}

/** Trim from both ends (copying)
 *
 * \param s string
 * \return trimmed string
 */
static inline std::string Ltrim_copy(std::string s) {
    Ltrim(s);
    return s;
}

/** Remove all tabs (in place)
 *
 * \param s string
 */
void Lsimplify(std::string & s)
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
bool LlisFloat(const string & str)
{
    std::istringstream iss(str);
    float f;
    iss >> noskipws >> f;
    return iss.eof() && !iss.fail();
}


MLookupManager * MLookupManager::lm = nullptr;

/** Returns instance of the Detector Manager class.
 *
 * \return manager instance
 */
MLookupManager * MLookupManager::instance()
{
    if (!lm)
        lm = new MLookupManager;

    return lm;
}

/** Shortcut
 * \return MLookupManager instance
 */
MLookupManager * lm()
{
    return MLookupManager::instance();
}

/** Default constructor
 */
MLookupManager::MLookupManager()
{
}

/** Destructor
 */
MLookupManager::~MLookupManager()
{
}

/** Parse source file
 *
 * \return success
 */
bool MLookupManager::parseSource()
{
    std::ifstream ifs(source);
    if (!ifs.is_open())
    {
        std::cerr << "Source " << source << " could not be open!" << std::endl;
        abort();
    }
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
        Ltrim(str);
        Lsimplify(str);

        size_t pos = 0;

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

                addLookupContainer(cont_name, new MLookupContainer(cont_name));

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
                containers[cont_name]->addLine(str);
                wn = WNContainerOrParam;
            }
        }
    }

    return true;
}

/** Print containers
 */
void MLookupManager::print() const
{
     std::map<std::string, MLookupContainer *>::const_iterator it = containers.begin();
     for (; it != containers.end(); ++it)
         it->second->print();
}

/** Add new parameter container.
 *
 * \param cont_name container name
 * \param parcont container object
 * \return success
 */
bool MLookupManager::addLookupContainer(const std::string& cont_name, MLookupContainer* lookupcont)
{
//     FIXME
//     MLookupContainer * lc = containers[cont_name];
//     if (!lc)
//     {
//         std::cerr << "Container " << cont_name << " doesn't exists!" << std::endl;
//         exit(EXIT_FAILURE);
//     }
// 
//     lookupcont->clear();
//     bool ret = lookupcont->getParams(lc);
//     if (!ret)
//     {
//         std::cerr << "Initialization of " << cont_name << " param container failed!" << std::endl;
//         exit(EXIT_FAILURE);
//     }
// 
    containers.insert(std::pair<std::string, MLookupContainer*>(cont_name, lookupcont));
    return true;
}

/** Get parameter container by name.
 *
 * \param cont_name container name
 * \return pointer to container
 */
MLookupContainer * MLookupManager::getLookupContainer(const std::string& cont_name)
{
    return containers[cont_name];
}
