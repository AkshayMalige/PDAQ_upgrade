// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdio>
#include <iostream>

#include "MLookupManager.h"
#include "MLookupContainer.h"
#include "MLookup.h"

/** \class MLookup
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

void MLookupChannel::setAddress(const char* address)
{
    sscanf(address, "%*s %*s %d %d %d\n", &m, &l, &s);
}

void MLookupChannel::print(const char * prefix)
{
    printf("%s %d %d %d\n", prefix, m, l, s);
}

MLookupBoard::MLookupBoard(UInt_t addr, UInt_t nchan) :
    addr(addr), nchan(nchan)
{
    channels = new MLookupChannel*[nchan];
    for (uint i = 0; i < nchan; ++i)
        channels[i] = 0;
}

MLookupBoard::~MLookupBoard()
{
    delete [] channels;
}

void MLookupBoard::print()
{
    char p[16];
    for (uint i = 0; i < nchan; ++i)
    {
        sprintf(p, "0x%x  %2d  ", addr, i);
        channels[i]->print(p);
    }
}


MLookupTable::MLookupTable(const std::string& container, UInt_t addr_min, UInt_t addr_max, UInt_t channels) :
    container(container), a_min(addr_min), a_max(addr_max), channels(channels), is_init(false)
{
    size_t nboards = addr_max - addr_min;
    boards = new MLookupBoard*[nboards];
    memset(boards, 0, sizeof(MLookupBoard*) * nboards);
}

void MLookupTable::init()
{
    MLookupContainer * lc = lm()->getLookupContainer(container);

    const LookupVector & lv = lc->getLines();

    size_t n = lv.size();
    for (uint i = 0; i < n; ++i)
    {
        uint addr, chan;
        char address[1024];
        sscanf(lv[i].c_str(), "%x %x %s", &addr, &chan, address);
        if (addr < a_min or addr > a_max)
        {
            std::cerr << "Can't add board " << addr << " inside (0x" << std::hex << a_min << ", 0x" << a_max << std::dec << "), skipping." << std::endl;
            continue;
        }
        uint idx = addr - a_min;
        if (boards[idx] == 0)
        {
            boards[idx] = new MLookupBoard(addr, channels);
            for (int i = 0; i < channels; ++i)
                boards[idx]->setChannel(i, initial());
        }

        boards[idx]->setAddress(chan, lv[i].c_str());
    }

    is_init = true;
}

MLookupTable::~MLookupTable()
{
    delete [] boards;
}

void MLookupTable::print()
{
    if (!is_init) init();

    printf("[%s]\n", container.c_str());
    size_t nboards = a_max - a_min;
    for (uint i = 0; i < nboards; ++i)
        if (boards[i]) boards[i]->print();
}
