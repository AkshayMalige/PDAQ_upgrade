// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef MLOOKUPCONTAINER_H
#define MLOOKUPCONTAINER_H

#include <TArrayI.h>
#include <TArrayF.h>
#include <TArrayD.h>

#include <vector>
#include <map>

//#include "MMAPTManager.h"

typedef std::vector<std::string> LookupVector;

class MLookupContainer
{
protected:
    // members
    std::string container;      ///< container name
    LookupVector lv;

public:
    // constructor
    MLookupContainer(const std::string & container);
    // destructor
    virtual ~MLookupContainer() {}

    void addLine(const std::string & line) { lv.push_back(line); }
    const LookupVector & getLines() const { return lv; }

    void print();
};

#endif // MLOOKUPCONTAINER_H
