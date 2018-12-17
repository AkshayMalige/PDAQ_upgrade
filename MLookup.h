// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef MLOOKUP_H
#define MLOOKUP_H

//#include "MMAPTManager.h"

class MLookupContainer;
class MLookupTable;

class MLookupChannel
{
public:
    UInt_t m, l, s, st;
    virtual void setAddress(const char * address);
    virtual void print(const char * prefix = 0);
};

class MLookupBoard
{
private:
    UInt_t addr, nchan;
    MLookupChannel ** channels;

public:
    MLookupBoard(UInt_t addr, UInt_t nchan);
    virtual ~MLookupBoard();

    void setAddress(UInt_t chan, const char * line) { channels[chan]->setAddress(line); }
    void setChannel(UInt_t chan, MLookupChannel * c) { channels[chan] = c; }
    MLookupChannel * getChannel(UInt_t chan) { return channels[chan]; }

    virtual void print();
};

class MLookupTable
{
protected:
    std::string container;
    int a_min, a_max, channels;
    bool is_init;

    MLookupBoard ** boards;

public:
    // constructor
    MLookupTable(const std::string & container, UInt_t addr_min, UInt_t addr_max, UInt_t channels = 49);
    // destructor
    virtual ~MLookupTable();

    virtual MLookupChannel * initial() = 0;

    MLookupChannel * getAddress(UInt_t addr, UInt_t chan) {
        if (!is_init) init();
	if (addr < a_min || addr > a_max) return nullptr;
        return boards[addr-a_min]->getChannel(chan);
    }

    virtual void print();

protected:
    void init();

};

#endif // MLOOKUP_H
