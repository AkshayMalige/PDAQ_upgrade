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

#ifndef MPARMANAGER_H
// @(#)lib/base:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#define MPARMANAGER_H

#include <map>

//#include "MMAPTManager.h"

class MParContainer;
class MPar;

class MParManager
{
protected:
    std::string source;             ///< Parameters source file
    std::string destination;        ///< Parameters destination file

    std::map<std::string, MParContainer *> containers;  ///< Containers
    std::map<std::string, MPar *> parconts;             ///< Parameters

    static MParManager * pm;        ///< Instance of the MParManager

private:
    // constructors
    MParManager();
    MParManager(MParManager const &) {}

    // methods
    /// Assignment operator
    /// \return this object
    MParManager & operator=(MParManager const &) { return *this; }

public:
    // instance method
    static MParManager * instance();
    // destructor
    ~MParManager();

    // methods
    /// Set parameters source
    /// \param source source file name
    void setParamSource(const std::string & source) { this->source = source; }
    /// Set parameters destination
    /// \param dest destination file name
    void setParamDest(const std::string & dest) { this->destination = dest; }

    bool parseSource();
    void writeDestination() const;

    bool addParameterContainer(const std::string & cont_name, MPar * parcont);
    MPar * getParameterContainer(const std::string & cont_name);

    void print() const;

private:
    /// Parser stepes
    enum WhatNext { WNContainer, WNContainerOrParam, WNParam, WNParamCont };
    WhatNext parseValues(const std::string & str, std::vector<std::string> & values);

};

extern MParManager * pm();

#endif // MPARMANAGER_H
