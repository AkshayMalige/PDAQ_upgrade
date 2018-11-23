// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include <sstream>

#include "MParContainer.h"

/** \class MParContainer
\ingroup lib_base

MPar is an abstract class to hold container and geometry parameters.

It must be derivated and pure virtual members defined.

The parameters are parsed from text file in MParManager and stored in the
MParContainer. The getParam() method reads content of the MParContainer and
fills variables inside the MPar object. The putParam method allows to update
parameters in the container and write to param file.

\sa MFibersStackCalibratorPar
\sa MFibersStackDigitizerPar
\sa MFibersStackGeomPar
*/

/** Constructor
 * \param container container name
 */
MParContainer::MParContainer(const std::string& container) : container(container)
{
}

/** Add key with integer value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, Int_t val)
{
    std::stringstream buff;
    buff << "  " << val;
    std::vector<std::string> v;
    v.push_back(buff.str());
    parameters[name] = TypeDataField("Int_t", v);

    return true;
}

/** Add key with float value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, Float_t val)
{
    std::stringstream buff;
    buff << "  " << val;
    std::vector<std::string> v;
    v.push_back(buff.str());
    parameters[name] = TypeDataField("Float_t", v);

    return true;
}

/** Add key with double precision float value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, Double_t val)
{
    std::stringstream buff;
    buff << "  " << val;
    std::vector<std::string> v;
    v.push_back(buff.str());
    parameters[name] = TypeDataField("Double_t", v);

    return true;
}

/** Add key with integer array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, const TArrayI & val)
{
    std::stringstream buff;
    std::vector<std::string> v;
    for (int i = 0; i < val.GetSize(); ++i)
    {
        buff << "  " << val[i];
        v.push_back(buff.str());
    }
    parameters[name] = TypeDataField("Int_t", v);

    return true;
}

/** Add key with float array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, const TArrayF & val)
{
   std::stringstream buff;
   std::vector<std::string> v;
    for (int i = 0; i < val.GetSize(); ++i)
    {
        buff << "  " << val[i];
        v.push_back(buff.str());
    }
    parameters[name] = TypeDataField("Float_t", v);

    return true;
}

/** Add key with double precision float array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::add(const std::string & name, const TArrayD & val)
{
    std::stringstream buff;
    std::vector<std::string> v;
    for (int i = 0; i < val.GetSize(); ++i)
    {
        buff << "  " << val[i];
        v.push_back(buff.str());
    }
    parameters[name] = TypeDataField("Double_t", v);

    return true;
}

/** Get key with integer value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, Int_t& val)
{
    ParMap::const_iterator it = parameters.find(name);
    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!!" << std::endl;
        return false;
    }
    if (it->second.first != "Int_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    std::stringstream buff(it->second.second[0]);
    buff >> val;

    return true;
}

/** Get key with float value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, Float_t& val)
{
    ParMap::const_iterator it = parameters.find(name);

    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!" << std::endl;
        return false;
    }
    if (it->second.first != "Float_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    std::stringstream buff(it->second.second[0]);
    buff >> val;

    return true;
}

/** Get key with double precision float value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, Double_t& val)
{
    ParMap::const_iterator it = parameters.find(name);

    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!" << std::endl;
        return false;
    }
    if (it->second.first != "Double_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    std::stringstream buff(it->second.second[0]);
    buff >> val;

    return true;
}

/** Get key with integer array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, TArrayI& val)
{
    ParMap::const_iterator it = parameters.find(name);

    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!" << std::endl;
        return false;
    }
    if (it->second.first != "Int_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    for (int i = 0; i < it->second.second.size(); ++i)
        val[i] = std::stoi(it->second.second[i].c_str());

    return true;
}

/** Get key with float array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, TArrayF& val)
{
    ParMap::const_iterator it = parameters.find(name);

    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!" << std::endl;
        return false;
    }
    if (it->second.first != "Float_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    for (int i = 0; i < it->second.second.size(); ++i)
        val[i] = std::stof(it->second.second[i].c_str());

    return true;
}

/** Get key with double precision float array value
 *
 * \param name key name
 * \param val value
 * \return success
 */
bool MParContainer::fill(const std::string & name, TArrayD& val)
{
    ParMap::const_iterator it = parameters.find(name);

    if (it == parameters.end())
    {
        std::cerr << "Parameter name " << name << " doesn't exists in the container!" << std::endl;
        return false;
    }
    if (it->second.first != "Double_t")
    {
        std::cerr << "Incorrect type for parameter name " << name << std::endl;
        return false;
    }

    for (int i = 0; i < it->second.second.size(); ++i)
        val[i] = std::stod(it->second.second[i].c_str());

    return true;
}

/** Print container
 */
void MParContainer::print()
{
    printf("Container [%s]\n", container.c_str());
    ParMap::const_iterator it = parameters.begin();
    for (; it != parameters.end(); ++it)
    {
        printf("%s:  %s", it->first.c_str(), it->second.first.c_str());
        for (int i = 0; i < it->second.second.size(); ++i)
        {
            printf("  %s", it->second.second[i].c_str());
        }
        putchar('\n');
    }
}

/** Init param with type and values
 *
 * \param name key name
 * \param type key type
 * \param values key values
 * \return success
 */
bool MParContainer::initParam(const std::string& name, const std::string& type, const std::vector<std::string> & values)
{
    parameters[name] = TypeDataField(type, values);
}
