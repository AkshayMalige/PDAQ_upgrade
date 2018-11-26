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
	sm_mods[i].fOffsetY = 0.0;
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
        //{cout<<"first layer "<<sm_mods[m].nLayers<<endl;
        return sm_mods[m].nLayers;//}
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
        //{cout<<"first offset z "<<sm_mods[m].fOffsetZ[l][p]<<endl;
        return sm_mods[m].fOffsetZ[l][p];//}
    else
        return -1;
}


Float_t MFTGeomPar::getOffsetY(Int_t m, Int_t l, Int_t p) const
{
    // return Z-coordinate of the beginning of the sublayer 's' in layer 'l'
    // m --module number
    // l -- layer number
    // s -- sublayer number
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
        //{cout<<"first offset z "<<sm_mods[m].fOffsetZ[l][p]<<endl;
        return sm_mods[m].fOffsetY[l][p];//}
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
    sm_mods[m].fOffsetY.ResizeTo(l, FWDET_STRAW_MAX_PLANES);

}

void MFTGeomPar::setStraws(Int_t m, Int_t s)
{
    // set number of straws for moduke 'm'
    // m -- module number
    // s -- number of straws
    if (m < getModules())
    {
        sm_mods[m].nStraws = s;
    }
}

void MFTGeomPar::setShortOffset(Int_t m, Int_t o)
{
    // set number of straws for moduke 'm'
    // m -- module number
    // o -- offset of short starw
    if (m < getModules())
    {
        sm_mods[m].nShortOffset = o;
    }
}

void MFTGeomPar::setShortWidth(Int_t m, Int_t w)
{
    // set number of straws for moduke 'm'
    // m -- module number
    // w -- width of short straws section
    if (m < getModules())
    {
        sm_mods[m].nShortWidth = w;
    }
}

void MFTGeomPar::setStrawRadius(Int_t m, Float_t r)
{
    // set straw radius in each sublayer of layer 'l'
    // l -- layer number
    // r -- straw radius
    if (m < getModules())
    {
        sm_mods[m].fStrawRadius = r;
    }
}

void MFTGeomPar::setStrawPitch(Int_t m, Float_t p)
{
    // set straws pitch in each sublayer of layer 'l'
    // l -- layer number
    // p -- straws pitch
    if (m < getModules())
    {
        sm_mods[m].fStrawPitch = p;
    }
}

void MFTGeomPar::setOffsetZ(Int_t m, Int_t l, Int_t p, Float_t z)
{
    // set Z-coordinate of the begining of 's' sublayer in layer 'l'
    // l -- layer number
    // s -- sublater number
    // z -- z-coordinate
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
    {
        sm_mods[m].fOffsetZ[l][p] = z;
    }
}

void MFTGeomPar::setOffsetX(Int_t m, Int_t l, Int_t p, Float_t x)
{
    // set X-coordinate of the begining of 's' sublayer in layer 'l'
    // l -- layer number
    // s -- sublater number
    // x -- z-coordinate
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
    {
        sm_mods[m].fOffsetX[l][p] = x;
    }
}

void MFTGeomPar::setOffsetY(Int_t m, Int_t l, Int_t p, Float_t y)
{
    // set Y-coordinate of the begining of 's' sublayer in layer 'l'
    // l -- layer number
    // s -- sublater number
    // y -- y-coordinate
    if (m < getModules() && l < getLayers(m) && p < FWDET_STRAW_MAX_PLANES)
    {
        sm_mods[m].fOffsetY[l][p] = y;
    }
}
void MFTGeomPar::setLayerRotation(Int_t m, Int_t l, Float_t r)
{
    // set lab transform for later 'l'
    // l -- layer number
    // lt -- lab transformation matrix [3x3 rot., 3-x transl.]
    if (m < getModules() && l < getLayers(m))
    {
        sm_mods[m].fLayerRotation[l] = r;
    }
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
    TArrayF par_offsetY(total_layers * FWDET_STRAW_MAX_PLANES);

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
		par_offsetY.SetAt(getOffsetY(i, l, s), 2*(cnt_layers+l) + s);

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
    parcont->add("fOffsetY",       par_offsetY);
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

    TArrayI par_straws(par_modules);
    if (!parcont->fill("nStraws", par_straws))
        return false;

    if (par_straws.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of straws does not fit to number of detectors");
        return false;
    }

    TArrayI par_shorto(par_modules);
    if (!parcont->fill("nShortOffset", par_shorto))
        return false;

    if (par_shorto.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of short straws offset does not fit to number of detectors");
        return false;
    }

    TArrayI par_shortw(par_modules);
    if (!parcont->fill("nShortWidth", par_shortw))
        return false;

    if (par_shortw.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of short straws section width does not fit to number of detectors");
        return false;
    }

    TArrayF par_strawRadius(par_modules);
    if (!parcont->fill("fStrawRadius", par_strawRadius))
        return false;

    if (par_strawRadius.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of strawRadius=%d does not fit to number of layers=%d", par_strawRadius.GetSize(), par_modules);
        return false;
    }

    TArrayF par_strawPitch(par_modules);
    if (!parcont->fill("fStrawPitch", par_strawPitch))
        return false;

    if (par_strawPitch.GetSize() != par_modules)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of strawPitch=%d does not fit to number of layers=%d", par_strawPitch.GetSize(), par_modules);
        return false;
    }

    Int_t cnt_layers = 0;

    for (Int_t d = 0; d < par_modules; ++d)
    {
        cnt_layers += par_layers[d];
    }
    const Int_t cnt_planes = cnt_layers * FWDET_STRAW_MAX_PLANES;

    TArrayF par_offsetZ(cnt_planes);
    if (!parcont->fill("fOffsetZ", par_offsetZ))
        return false;

    if (par_offsetZ.GetSize() != cnt_planes)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of planeZ=%d does not fit to number of planes=%d", par_offsetZ.GetSize(), cnt_planes);
        return false;
    }

    TArrayF par_offsetX(cnt_planes);
    if (!parcont->fill("fOffsetX", par_offsetX))
        return false;

    if (par_offsetX.GetSize() != cnt_planes)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of planeX=%d does not fit to number of planes=%d", par_offsetX.GetSize(), cnt_planes);
        return false;
    }
    
    TArrayF par_offsetY(cnt_planes);
    if (!parcont->fill("fOffsetY", par_offsetY))
        return false;

    if (par_offsetY.GetSize() != cnt_planes)
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of planeY=%d does not fit to number of planes=%d", par_offsetY.GetSize(), cnt_planes);
        return false;
    }

    TArrayF par_layerRotation(total_layers);
    if (!parcont->fill("fLayerRotation", par_layerRotation))
        return false;

    if (par_layerRotation.GetSize() != (Int_t)(total_layers))
    {
        Error("HFwDetStrawGeomPar::getParams(HParamList* l)",
              "Array size of layerRotation=%d does not fit to number of layers=%d", par_layerRotation.GetSize(), total_layers);
        return false;
    }

     cnt_layers = 0;

     setModules(par_modules);

    for (Int_t i = 0; i < par_modules; ++i)
    {
    //     // get number of layers
        Int_t layers = par_layers[i];

        // set number of layers
        setLayers(i, layers);
        setStrawRadius(i, par_strawRadius[i]);
        setStrawPitch(i, par_strawPitch[i]);

        // iterate over layers
        for (Int_t l = 0; l < layers; ++l)
        {
            for (Int_t s = 0; s < FWDET_STRAW_MAX_PLANES; ++s)
            {
                setOffsetZ(i, l, s, par_offsetZ[2*(cnt_layers + l) + s]);
                setOffsetX(i, l, s, par_offsetX[2*(cnt_layers + l) + s]);
		setOffsetY(i, l, s, par_offsetY[2*(cnt_layers + l) + s]);

            }
            setLayerRotation(i, l, par_layerRotation[cnt_layers + l]);
        }

        // set number of straws in each block
        setStraws(i, par_straws[i]);
        // set short straws properties
        setShortOffset(i, par_shorto[i]);
        setShortWidth(i, par_shortw[i]);

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
    for (int m = 0; m < nModules; ++m)
    {
        printf(" +++\n layers = %d\n", sm_mods[m].nLayers);
//        printf(" +++\n layers = %d\n", mods[m].layers);
        printf(" Straws:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2d", sm_mods[m].nStraws);

        printf("\n layrot:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2.f", sm_mods[m].fLayerRotation[l]);

        printf("\n  off x:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            {for (int p=0; p< FWDET_STRAW_MAX_PLANES; ++p)
            {printf(" %2.3f", sm_mods[m].fOffsetX[l][p]);}}

        printf("\n  off Z:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            {for (int p =0; p< FWDET_STRAW_MAX_PLANES; ++p)
                {printf(" %2.3f", sm_mods[m].fOffsetZ[l][p]);}}
                
        printf("\n  off Y:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            {for (int p =0; p< FWDET_STRAW_MAX_PLANES; ++p)
                {printf(" %2.3f", sm_mods[m].fOffsetY[l][p]);}}

        printf("\n  short offset:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2d", sm_mods[m].nShortOffset);
        printf("\n  short width:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2d", sm_mods[m].nShortWidth);

        printf("\n  radius:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2.3f", sm_mods[m].fStrawRadius);
        printf("\n  pitch:");
        for (int l = 0; l < sm_mods[m].nLayers; ++l)
            printf(" %2.3f", sm_mods[m].fStrawPitch);
        putchar('\n');
    }
}


