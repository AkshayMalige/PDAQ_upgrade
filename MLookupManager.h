/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  <copyright holder> <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef MLOOKUPMANAGER_H
// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#define MLOOKUPMANAGER_H

#include <map>

//#include "MMAPTManager.h"

class MLookupContainer;
class MLookup;

class MLookupManager
{
protected:
    std::string source;             ///< Parameters source file
    std::string destination;        ///< Parameters destination file

    std::map<std::string, MLookupContainer *> containers;  ///< Containers

    static MLookupManager * lm;        ///< Instance of the MLookupManager

private:
    // constructors
    MLookupManager();
    MLookupManager(MLookupManager const &) {}

    // methods
    /// Assignment operator
    /// \return this object
    MLookupManager & operator=(MLookupManager const &) { return *this; }

public:
    // instance method
    static MLookupManager * instance();
    // destructor
    ~MLookupManager();

    // methods
    /// Set parameters source
    /// \param source source file name
    void setSource(const std::string & source) { this->source = source; }

    bool parseSource();

    bool addLookupContainer(const std::string & cont_name, MLookupContainer * lookupcont);
    MLookupContainer * getLookupContainer(const std::string & cont_name);

    void print() const;

private:
    /// Parser stepes
    enum WhatNext { WNContainer, WNContainerOrParam, WNParam, WNParamCont };

};

extern MLookupManager * lm();

#endif // MLOOKUPMANAGER_H
