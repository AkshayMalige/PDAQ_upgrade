// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#include "MPar.h"

#include "MParContainer.h"

/** \class MPar
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

/** Default constructor
 */
MPar::MPar() : parcont(nullptr)
{
}

/** Destructor
 */
MPar::~MPar()
{
}

/** Print container info
 */
void MPar::print() const
{
    parcont->print();
}
