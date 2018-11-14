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

#include "FTGeo.h"


#include "MParContainer.h"

/** \class MFTGeomPar
\ingroup lib_fibers_stack

A container for Fibers Stack geometry parameters

\sa MPar

*/
using namespace std;
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
    nModules = 0;

    for (Int_t i = 0; i < FWDET_STRAW_MAX_MODULES; ++i)
    {
        // clear all variables
        sm_mods[i].nLayers = 0;
        sm_mods[i].nStraws = 0;
        sm_mods[i].nShortOffset = 0;
        sm_mods[i].nShortWidth = 0;
        sm_mods[i].fStrawRadius = 0.0;
        sm_mods[i].fStrawPitch = 0.0;
        sm_mods[i].fOffsetZ = 0.0;
        sm_mods[i].fOffsetX = 0.0;
        sm_mods[i].fLayerRotation = 0;
    }
}

// Int_t MFTGeomPar::getModules() const
// {
//     // return number of layers in single detector
//     return nModules;
// }

/** Get parameters
 * \sa MPar::getParams()
 * \param parcont pointer to container object
 * \return success
 */

Int_t MFTGeomPar::getModules() const
{
    // return number of layers in single detector
    return nModules;
}

Int_t MFTGeomPar::getLayers(Int_t m) const
{
    // return number of layers in module 'm'
    // m -- module number
    if (m < getModules())
        return sm_mods[m].nLayers;
    else
        return -1;
}

Int_t MFTGeomPar::getStraws(Int_t m) const
{
    // return number of straws in single layer of module 'm'
    // l -- layer number
    if (m < getModules())
        return sm_mods[m].nStraws;//[_off_straws[l][d][b]];
    else
        return -1;
}

Int_t MFTGeomPar::getShortWidth(Int_t m) const
{
    // return number of straws in single layer of module 'm'
    // l -- layer number
    if (m < getModules())
        return sm_mods[m].nShortWidth;//[_off_straws[l][d][b]];
    else
        return -1;
}

Float_t MFTGeomPar::getLayerRotation(Int_t m, Int_t l) const
{
    // return transformation matrix for a layer 'l'
    // l -- layer number
    if (m < getModules() && l < getLayers(m))
        return sm_mods[m].fLayerRotation[l];
    else
        return -1;
}


Int_t MFTGeomPar::getShortOffset(Int_t m) const
{
    // return number of straws in single layer of module 'm'
    // l -- layer number
    if (m < getModules())
        return sm_mods[m].nShortOffset;//[_off_straws[l][d][b]];
    else
        return -1;
}


Float_t MFTGeomPar::getOffsetX(Int_t m, Int_t l, Int_t p) const
{
    // return X-coordinate of the beginning of the sublayer 's' in layer 'l'
    // m -- module number
    // l -- layer number
    // s -- sublayer number
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
        return sm_mods[m].fOffsetX[l][p];
    else
        return -1;
}


Float_t MFTGeomPar::getStrawRadius(Int_t m) const
{
    // return straw radius in single sublayer of layer 'l'
    // l -- layer number
    if (m < getModules())
        return sm_mods[m].fStrawRadius;
    else
        return -1;
}

Float_t MFTGeomPar::getStrawPitch(Int_t m) const
{
    // return straw pitch in single sublayer of layer 'l'
    // l -- layer number
    if (m < getModules())
        return sm_mods[m].fStrawPitch;
    else
        return -1;
}

Float_t MFTGeomPar::getOffsetZ(Int_t m, Int_t l, Int_t p) const
{
    // return Z-coordinate of the beginning of the sublayer 's' in layer 'l'
    // m --module number
    // l -- layer number
    // s -- sublayer number
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
        return sm_mods[m].fOffsetZ[l][p];
    else
        return -1;
}


void MFTGeomPar::setModules(int m)
{
    if (m <= FWDET_STRAW_MAX_MODULES)
        nModules = m;
}


void MFTGeomPar::setLayers(Int_t m, Int_t l)
{
    // set number of layers, this function automatically
    // resizes all depending arrays
    sm_mods[m].nLayers = l;
    sm_mods[m].fLayerRotation.Set(l);
    sm_mods[m].fOffsetZ.ResizeTo(l, FWDET_STRAW_MAX_PLANES);
    sm_mods[m].fOffsetX.ResizeTo(l, FWDET_STRAW_MAX_PLANES);
}

bool MFTGeomPar::putParams(MParContainer* parcont) const
{
     if (!parcont) return 0;

      Int_t total_layers = 0;

    // first find total number of layers
    for (Int_t i = 0; i < nModules; ++i)
    {
        total_layers += getLayers(i);
    }  

    TArrayI par_layers(nModules);
    TArrayI par_straws(nModules);
    TArrayI par_shorto(nModules);
    TArrayI par_shortw(nModules);

    TArrayF par_offsetX(total_layers * FWDET_STRAW_MAX_PLANES);
    TArrayF par_offsetZ(total_layers * FWDET_STRAW_MAX_PLANES);

    TArrayF par_strawR(nModules);
    TArrayF par_strawP(nModules);

    TArrayF par_layerRotation(total_layers);

    Int_t cnt_layers = 0; 

    for (Int_t i = 0; i < nModules; ++i)
    {
        // get number of layers
        Int_t layers = getLayers(i);

        // set number of layers
        par_layers.SetAt(layers, i);

        par_strawR.SetAt(getStrawRadius(i), i);
        par_strawP.SetAt(getStrawPitch(i), i);

        // iterate over layers
        for (Int_t l = 0; l < layers; ++l)
        {
            for (Int_t s = 0; s < FWDET_STRAW_MAX_PLANES; ++s)
            {
                par_offsetZ.SetAt(getOffsetZ(i, l, s), 2*(cnt_layers+l) + s);
                par_offsetX.SetAt(getOffsetX(i, l, s), 2*(cnt_layers+l) + s);
            }

            par_layerRotation.SetAt(getLayerRotation(i, l), cnt_layers + l);
        }

        // set number of straws in each block
        par_straws.SetAt(getStraws(i), i);
        // set short straws values
        par_shorto.SetAt(getShortOffset(i), i);
        par_shortw.SetAt(getShortWidth(i), i);

        cnt_layers += layers;
    }

    parcont->add("nModules",       nModules);
    parcont->add("nLayers",        par_layers);
    parcont->add("nStraws",        par_straws);
    parcont->add("nShortOffset",   par_shorto);
    parcont->add("nShortWidth",    par_shortw);
    parcont->add("fOffsetZ",       par_offsetZ);
    parcont->add("fOffsetX",       par_offsetX);
    parcont->add("fStrawRadius",   par_strawR);
    parcont->add("fStrawPitch",    par_strawP);
    parcont->add("fLayerRotation", par_layerRotation);

}


bool MFTGeomPar::getParams(MParContainer* parcont)
{
// gets the parameters from the list (read from input)
    if (!parcont) return false;

    Int_t par_modules;
    if (!parcont->fill("nModules", par_modules))

        return false;


    TArrayI par_layers(par_modules);
    if (!parcont->fill("nLayers", par_layers))
        return false;

    if (par_layers.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(MFTGeomPar* parcont)",
              "Array size of layers does not fit to number of detectors");
        return false;
    }


    Int_t total_layers = 0;
    
    for (Int_t d = 0; d < par_modules; ++d)
    {
        total_layers += par_layers[d];
    }

    // TArrayI par_straws;
    // if (!parcont->fill("nStraws", par_straws))
    //     return false;

    // if (par_straws.GetSize() != par_modules)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of straws does not fit to number of detectors");
    //     return false;
    // }

    // TArrayI par_shorto;
    // if (!parcont->fill("nShortOffset", &par_shorto))
    //     return false;

    // if (par_shorto.GetSize() != par_modules)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of short straws offset does not fit to number of detectors");
    //     return false;
    // }

    // TArrayI par_shortw;
    // if (!parcont->fill("nShortWidth", &par_shortw))
    //     return false;

    // if (par_shortw.GetSize() != par_modules)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of short straws section width does not fit to number of detectors");
    //     return false;
    // }

    // TArrayF par_strawRadius;
    // if (!parcont->fill("fStrawRadius", &par_strawRadius))
    //     return false;

    // if (par_strawRadius.GetSize() != par_modules)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of strawRadius=%d does not fit to number of layers=%d", par_strawRadius.GetSize(), par_modules);
    //     return false;
    // }

    // TArrayF par_strawPitch;
    // if (!parcont->fill("fStrawPitch", &par_strawPitch))
    //     return false;

    // if (par_strawPitch.GetSize() != par_modules)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of strawPitch=%d does not fit to number of layers=%d", par_strawPitch.GetSize(), par_modules);
    //     return false;
    // }

     Int_t cnt_layers = 0;
    // for (Int_t d = 0; d < par_modules; ++d)
    // {
    //     cnt_layers += par_layers[d];
    // }
    // const Int_t cnt_planes = cnt_layers * FWDET_STRAW_MAX_PLANES;

    // TArrayF par_offsetZ(cnt_planes);
    // if (!parcont->fill("fOffsetZ", &par_offsetZ))
    //     return false;

    // if (par_offsetZ.GetSize() != cnt_planes)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of planeZ=%d does not fit to number of planes=%d", par_offsetZ.GetSize(), cnt_planes);
    //     return false;
    // }

    // TArrayF par_offsetX(cnt_planes);
    // if (!parcont->fill("fOffsetX", &par_offsetX))
    //     return false;

    // if (par_offsetX.GetSize() != cnt_planes)
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of planeX=%d does not fit to number of planes=%d", par_offsetX.GetSize(), cnt_planes);
    //     return false;
    // }

    // TArrayF par_layerRotation(total_layers);
    // if (!parcont->fill("fLayerRotation", &par_layerRotation))
    //     return false;

    // if (par_layerRotation.GetSize() != (Int_t)(total_layers))
    // {
    //     Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
    //           "Array size of layerRotation=%d does not fit to number of layers=%d", par_layerRotation.GetSize(), total_layers);
    //     return false;
    // }

     cnt_layers = 0;

     setModules(par_modules);

    for (Int_t i = 0; i < par_modules; ++i)
    {
    //     // get number of layers
        Int_t layers = par_layers[i];

        // set number of layers
        setLayers(i, layers);

    //     setStrawRadius(i, par_strawRadius[i]);
    //     setStrawPitch(i, par_strawPitch[i]);

    //     // iterate over layers
    //     for (Int_t l = 0; l < layers; ++l)
    //     {
    //         for (Int_t s = 0; s < FWDET_STRAW_MAX_PLANES; ++s)
    //         {
    //             setOffsetZ(i, l, s, par_offsetZ[2*(cnt_layers + l) + s]);
    //             setOffsetX(i, l, s, par_offsetX[2*(cnt_layers + l) + s]);
    //         }
    //         setLayerRotation(i, l, par_layerRotation[cnt_layers + l]);
    //     }

    //     // set number of straws in each block
    //     setStraws(i, par_straws[i]);
    //     // set short straws properties
    //     setShortOffset(i, par_shorto[i]);
    //     setShortWidth(i, par_shortw[i]);

        cnt_layers += layers;
    }

    return true;
}

/** Put parameters
 * \sa MPar::putParams()
 * \param parcont pointer to container object
 * \return success
 */


/** Print parameters
 */
void MFTGeomPar::print() const
{
    printf("Number of modules = %d\n", nModules);
    // for (int m = 0; m < modules; ++m)
    // {
    //     printf(" +++\n layers = %d\n", mods[m].layers);
    //     printf(" fibers:");
    //     for (int l = 0; l < mods[m].layers; ++l)
    //         printf(" %2d", mods[m].fibers[l]);
    //     printf("\n layrot:");
    //     for (int l = 0; l < mods[m].layers; ++l)
    //         printf(" %2.0f", mods[m].layer_rotation[l]);
    //     printf("\n  off x:");
    //     for (int l = 0; l < mods[m].layers; ++l)
    //         printf(" %2.0f", mods[m].fiber_offset_x[l]);
    //     printf("\n  off y:");
    //     for (int l = 0; l < mods[m].layers; ++l)
    //         printf(" %2.0f", mods[m].fiber_offset_y[l]);
    //     printf("\n  pitch:");
    //     for (int l = 0; l < mods[m].layers; ++l)
    //         printf(" %2.0f", mods[m].fibers_pitch[l]);
    //     putchar('\n');
    // }
}

// /** Get number of layers in the module
//  * \param m module
//  * \return number of layers
//  */
// Int_t MFTGeomPar::getLayers(Int_t m) const
// {
//     if (mods and m < modules)
//         return mods[m].layers;
//     else
//         return -1;
// }

// /** Get number of fibers in the layer
//  * \param m module
//  * \param l layer
//  * \return number of fibers
//  */
// Int_t MFTGeomPar::getFibers(Int_t m, Int_t l) const
// {
//     if (mods and m < modules and l < mods[m].layers)
//         return mods[m].fibers[l];
//     else
//         return -1;
// }

// /** Get layer rotation
//  * \param m module
//  * \param l layer
//  * \return layer rotation
//  */
// Float_t MFTGeomPar::getLayerRotation(Int_t m, Int_t l) const
// {
//     if (mods and m < modules and l < mods[m].layers)
//         return mods[m].layer_rotation[l];
//     else
//         return -10000.;
// }

// /** Get fibers offset X
//  * \param m module
//  * \param l layer
//  * \return offset X
//  */
// Float_t MFTGeomPar::getFiberOffsetX(Int_t m, Int_t l) const
// {
//     if (mods and m < modules and l < mods[m].layers)
//         return mods[m].fiber_offset_x[l];
//     else
//         return -10000.;
// }

// /** Get layers offset Y
//  * \param m module
//  * \param l layer
//  * \return offset Y
//  */
// Float_t MFTGeomPar::getFiberOffsetY(Int_t m, Int_t l) const
// {
//     if (mods and m < modules and l < mods[m].layers)
//         return mods[m].fiber_offset_y[l];
//     else
//         return -10000.;
// }

// /** Get fibers pitch in a layer
//  * \param m module
//  * \param l layer
//  * \return fibers pitch
//  */
// Float_t MFTGeomPar::getFibersPitch(Int_t m, Int_t l) const
// {
//     if (mods and m < modules and l < mods[m].layers)
//         return mods[m].fibers_pitch[l];
//     else
//         return -10000.;
// }
