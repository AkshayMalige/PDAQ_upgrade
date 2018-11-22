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

#include "MLookupContainer.h"

/** \class MLookupContainer
\ingroup lib_base

MLookup is an abstract class to hold container and geometry parameters.

It must be derivated and pure virtual members defined.

The parameters are parsed from text file in MLookupManager and stored in the
MLookupContainer. The getParam() method reads content of the MLookupContainer and
fills variables inside the MLookup object. The putParam method allows to update
parameters in the container and write to param file.

\sa MFibersStackCalibratorPar
\sa MFibersStackDigitizerPar
\sa MFibersStackGeomPar
*/

/** Constructor
 * \param container container name
 */
MLookupContainer::MLookupContainer(const std::string& container)
    : container(container)
{
}

/** Print container
 */
void MLookupContainer::print()
{
    printf("Lookup Table [%s]\n", container.c_str());
    LookupVector::const_iterator it = lv.begin();
    for (; it != lv.end(); ++it)
    {
        printf("%s", it->c_str());
    }
}
