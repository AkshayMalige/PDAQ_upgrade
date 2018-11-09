// @(#)lib/fibers_stack:$Id$
// Author: Rafal Lalik  18/11/2017

/*************************************************************************
 * Copyright (C) 2017-2018, Rafał Lalik.                                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $MAPTSYS/LICENSE.                         *
 * For the list of contributors see $MAPTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>

#include "MFTGeomPar.h"

#include "MParContainer.h"

/** \class MFTGeomPar
\ingroup lib_fibers_stack

A container for Fibers Stack geometry parameters

\sa MPar

*/

/** Constructor
 */
MFTGeomPar::MFTGeomPar() : MPar(), mods(nullptr)
{
}

/** Destructor
 */
MFTGeomPar::~MFTGeomPar()
{
    delete [] mods;
}

/** Clear parameters
 */
void MFTGeomPar::clear()
{
    delete [] mods;
    mods = nullptr;
    modules = 0;
}

/** Get parameters
 * \sa MPar::getParams()
 * \param parcont pointer to container object
 * \return success
 */
bool MFTGeomPar::getParams(MParContainer* parcont)
{
    if (!parcont->fill("nModules", modules)) return false;

    if (mods) delete [] mods;
    mods = new SingleModule[modules];

    // get layers
    TArrayI _l(modules);
    if (!parcont->fill("nLayers", _l)) return false;
    if (_l.GetSize() != modules)
    {
        std::cerr << "Size of nLayers doesn't match nModules" << std::endl;
        return false;
    }

    int n_layers = _l.GetSum();

    // get fibers
    TArrayI _f(modules*n_layers);
    if (!parcont->fill("nFibers", _f)) return false;
    if (_f.GetSize() != (modules*n_layers))
    {
        std::cerr << "Size of nFibers doesn't match nModules*nLayers" << std::endl;
        return false;
    }

    TArrayF _lr(modules*n_layers);
    if (!parcont->fill("fLayerRotation", _lr)) return false;
    if (_lr.GetSize() != (modules*n_layers))
    {
        std::cerr << "Size of fLayerRotation doesn't match nModules" << std::endl;
        return false;
    }

    TArrayF _fox(modules*n_layers);
    if (!parcont->fill("fFiberOffsetX", _fox)) return false;
    if (_fox.GetSize() != (modules*n_layers))
    {
        std::cerr << "Size of fFiberOffsetX doesn't match nModules" << std::endl;
        return false;
    }

    TArrayF _foy(modules*n_layers);
    if (!parcont->fill("fFiberOffsetY", _foy)) return false;
    if (_foy.GetSize() != (modules*n_layers))
    {
        std::cerr << "Size of fFiberOffsetY doesn't match nModules" << std::endl;
        return false;
    }

    TArrayF _fsp(modules*n_layers);
    if (!parcont->fill("fFibersPitch", _fsp)) return false;
    if (_fsp.GetSize() != (modules*n_layers))
    {
        std::cerr << "Size of fFibersPitch doesn't match nModules" << std::endl;
        return false;
    }

    for (int m = 0; m < modules; ++m)
    {
        mods[m].layers = _l[m];
        mods[m].fibers.Set(mods[m].layers);
        mods[m].layer_rotation.Set(mods[m].layers);
        mods[m].fiber_offset_x.Set(mods[m].layers);
        mods[m].fiber_offset_y.Set(mods[m].layers);
        mods[m].fibers_pitch.Set(mods[m].layers);

        for (int l = 0; l < _l[m]; ++l)
        {
            Int_t fibers = _f[n_layers * m + l];
            Int_t layrot = _lr[n_layers * m + l];
            Float_t fox = _fox[n_layers * m + l];
            Float_t foy = _foy[n_layers * m + l];
            Float_t fsp = _fsp[n_layers * m + l];

            mods[m].fibers[l] = fibers;
            mods[m].layer_rotation[l] = layrot;
            mods[m].fiber_offset_x[l] = fox;
            mods[m].fiber_offset_y[l] = foy;
            mods[m].fibers_pitch[l] = fsp;
        }
    }

    return true;
}

/** Put parameters
 * \sa MPar::putParams()
 * \param parcont pointer to container object
 * \return success
 */
bool MFTGeomPar::putParams(MParContainer* parcont) const
{
}

/** Print parameters
 */
void MFTGeomPar::print() const
{
    printf("Number of modules = %d\n", modules);
    for (int m = 0; m < modules; ++m)
    {
        printf(" +++\n layers = %d\n", mods[m].layers);
        printf(" fibers:");
        for (int l = 0; l < mods[m].layers; ++l)
            printf(" %2d", mods[m].fibers[l]);
        printf("\n layrot:");
        for (int l = 0; l < mods[m].layers; ++l)
            printf(" %2.0f", mods[m].layer_rotation[l]);
        printf("\n  off x:");
        for (int l = 0; l < mods[m].layers; ++l)
            printf(" %2.0f", mods[m].fiber_offset_x[l]);
        printf("\n  off y:");
        for (int l = 0; l < mods[m].layers; ++l)
            printf(" %2.0f", mods[m].fiber_offset_y[l]);
        printf("\n  pitch:");
        for (int l = 0; l < mods[m].layers; ++l)
            printf(" %2.0f", mods[m].fibers_pitch[l]);
        putchar('\n');
    }
}

/** Get number of layers in the module
 * \param m module
 * \return number of layers
 */
Int_t MFTGeomPar::getLayers(Int_t m) const
{
    if (mods and m < modules)
        return mods[m].layers;
    else
        return -1;
}

/** Get number of fibers in the layer
 * \param m module
 * \param l layer
 * \return number of fibers
 */
Int_t MFTGeomPar::getFibers(Int_t m, Int_t l) const
{
    if (mods and m < modules and l < mods[m].layers)
        return mods[m].fibers[l];
    else
        return -1;
}

/** Get layer rotation
 * \param m module
 * \param l layer
 * \return layer rotation
 */
Float_t MFTGeomPar::getLayerRotation(Int_t m, Int_t l) const
{
    if (mods and m < modules and l < mods[m].layers)
        return mods[m].layer_rotation[l];
    else
        return -10000.;
}

/** Get fibers offset X
 * \param m module
 * \param l layer
 * \return offset X
 */
Float_t MFTGeomPar::getFiberOffsetX(Int_t m, Int_t l) const
{
    if (mods and m < modules and l < mods[m].layers)
        return mods[m].fiber_offset_x[l];
    else
        return -10000.;
}

/** Get layers offset Y
 * \param m module
 * \param l layer
 * \return offset Y
 */
Float_t MFTGeomPar::getFiberOffsetY(Int_t m, Int_t l) const
{
    if (mods and m < modules and l < mods[m].layers)
        return mods[m].fiber_offset_y[l];
    else
        return -10000.;
}

/** Get fibers pitch in a layer
 * \param m module
 * \param l layer
 * \return fibers pitch
 */
Float_t MFTGeomPar::getFibersPitch(Int_t m, Int_t l) const
{
    if (mods and m < modules and l < mods[m].layers)
        return mods[m].fibers_pitch[l];
    else
        return -10000.;
}
