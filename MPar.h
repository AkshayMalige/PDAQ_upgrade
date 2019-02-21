// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef MPAR_H
#define MPAR_H

//#include "MMAPTManager.h"

class MParContainer;

class MPar
{
    friend class MParManager; ///< MParManager is friend
  protected:
    // members
    MParContainer* parcont; ///< Pointer to object of parameter container

  public:
    // constructor
    MPar();
    // destructor
    virtual ~MPar();

    // methods
    /// Clear parameters
    virtual void clear() = 0;
    virtual void print() const;

  protected:
    /// Get parameters from container
    /// \param parcont pointer to container object
    /// \return success
    virtual bool getParams(MParContainer* parcont) = 0;
    /// Put parameters into container
    /// \param parcont pointer to container object
    /// \return success
    virtual bool putParams(MParContainer* parcont) const = 0;
};

#endif // MPAR_H
