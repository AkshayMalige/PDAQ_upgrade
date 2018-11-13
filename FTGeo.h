// @(#)lib/fibers_stack:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef MFTGEOMPAR_H
#define MFTGEOMPAR_H

#include "FTdef.h"
#include <TArrayI.h>
#include <TArrayF.h>
#include "TMatrixF.h"
#include "MPar.h"

class MFTGeomPar : public MPar
{
protected:
    // members
    Int_t nModules;                  ///< number of modules
    struct SingleModule             ///< single module configuration
    {
        Int_t nLayers;               ///< numbre of layers
        TArrayI nStraws; 
        Int_t    nShortOffset;      // offset of a short straw
        Int_t    nShortWidth; 

        Float_t  fStrawRadius;   // [layer]
        Float_t  fStrawPitch;    // [layer]

        TMatrixF fOffsetZ;     ///< offset of the first fiber in the layer
        TMatrixF fOffsetX;     ///< offset of the layers

        TArrayF fLayerRotation;     ///< layer rotation around the axis
        TArrayF fibers_pitch;       ///< fibers pitch in a layer
    }sm_mods[FWDET_STRAW_MAX_MODULES];

    SingleModule * mods;            ///< params for each module

public:
    MFTGeomPar();
    virtual ~MFTGeomPar();

    bool getParams(MParContainer * parcont);
    bool putParams(MParContainer * parcont) const;
    void clear();
    void print() const;

    /// Get number of modules
    /// \return number of modules
    Int_t getModules() const { return nModules; }
    Int_t getLayers(Int_t m) const;
    Int_t getStraws(Int_t m, Int_t l) const;
    Int_t getShortOffset(Int_t m) const;
    Int_t getShortWidth(Int_t m) const;
    Float_t getStrawRadius(Int_t m) const;
    Float_t getStrawPitch(Int_t m) const;
    Float_t getOffsetZ(Int_t m, Int_t l, Int_t p) const;
    Float_t getOffsetX(Int_t m, Int_t l, Int_t p) const;
    Float_t getLayerRotation(Int_t m, Int_t l) const;

    void setModules(Int_t m);
    void setLayers(Int_t m, Int_t l);
    void setStraws(Int_t m, Int_t s);
    void setShortOffset(Int_t m, Int_t o);
    void setShortWidth(Int_t m, Int_t w);
    void setStrawRadius(Int_t m, Float_t r);
    void setStrawPitch(Int_t m, Float_t p);
    void setOffsetZ(Int_t m, Int_t l, Int_t p, Float_t z);
    void setOffsetX(Int_t m, Int_t l, Int_t p, Float_t x);
    void setLayerRotation(Int_t m, Int_t l, Float_t r);
    
};

#endif // MFIBERSSTACKGEOMPAR_H
