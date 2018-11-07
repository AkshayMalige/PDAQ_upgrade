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
