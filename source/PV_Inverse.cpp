/*
 SuperCollider real time audio synthesis system
 Copyright (c) 2002 James McCartney. All rights reserved.
 http://www.audiosynth.com
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */



#include "FFT_UGens.h"


InterfaceTable *ft;

struct PV_OutOfPlace : Unit
{
    int m_numbins;
    float *m_tempbuf;
};

struct PV_DistrN : PV_OutOfPlace
{
};

struct PV_Inverse_0 : PV_OutOfPlace
{
};

struct PV_Inverse_1 : PV_OutOfPlace
{
};

struct PV_Inverse_2 : PV_OutOfPlace
{
};

struct PV_Inverse_Exp : PV_OutOfPlace
{
};

struct PV_Inv : PV_OutOfPlace
{
};

struct PV_Transp : PV_OutOfPlace
{
};

struct PV_InvTransp : PV_OutOfPlace
{
};



//////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
{
    void PV_DistrN_Ctor(PV_DistrN *unit);
    void PV_DistrN_Dtor(PV_DistrN *unit);
    void PV_DistrN_next(PV_DistrN *unit, int inNumSamples);
    
    void PV_Inverse_0_Ctor(PV_Inverse_0 *unit);
    void PV_Inverse_0_Dtor(PV_Inverse_0 *unit);
    void PV_Inverse_0_next(PV_Inverse_0 *unit, int inNumSamples);
    
    void PV_Inverse_1_Ctor(PV_Inverse_1 *unit);
    void PV_Inverse_1_Dtor(PV_Inverse_1 *unit);
    void PV_Inverse_1_next(PV_Inverse_1 *unit, int inNumSamples);
    
    void PV_Inverse_2_Ctor(PV_Inverse_2 *unit);
    void PV_Inverse_2_Dtor(PV_Inverse_2 *unit);
    void PV_Inverse_2_next(PV_Inverse_2 *unit, int inNumSamples);
    
    void PV_Inverse_Exp_Ctor(PV_Inverse_Exp *unit);
    void PV_Inverse_Exp_Dtor(PV_Inverse_Exp *unit);
    void PV_Inverse_Exp_next(PV_Inverse_Exp *unit, int inNumSamples);
    
    void PV_Inv_Ctor(PV_Inv *unit);
    void PV_Inv_Dtor(PV_Inv *unit);
    void PV_Inv_next(PV_Inv *unit, int inNumSamples);
    
    void PV_Transp_Ctor(PV_Transp *unit);
    void PV_Transp_Dtor(PV_Transp *unit);
    void PV_Transp_next(PV_Transp *unit, int inNumSamples);
    
    void PV_InvTransp_Ctor(PV_InvTransp *unit);
    void PV_InvTransp_Dtor(PV_InvTransp *unit);
    void PV_InvTransp_next(PV_InvTransp *unit, int inNumSamples);
}


///////////////////////////////////////////////////////////////////////////////////////////////

//***********************
//PV_DistrN
//***********************

void PV_DistrN_next(PV_DistrN *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        m_k = (float)(i + 1.f);
        
        //ajustements
        m_k = std::floor( numbins / m_k);
        m_k = std::max( m_k , 1.f );
        m_k = std::floor( numbins / m_k);
        m_k = std::max( m_k , 1.f );
        
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_DistrN_Ctor(PV_DistrN *unit)
{
    SETCALC(PV_DistrN_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_DistrN_Dtor(PV_DistrN *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}


///////////////////////////////////////////////////////////////////////////////////////////////
//*************
//PV_Inverse_0 — algorithme original (round - not safe / glitch !)
//*************

void PV_Inverse_0_next(PV_Inverse_0 *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr, freqC, freqT,fond0;
    float krC,krT;
    int size;
    
    freqC = ZIN0(1);
    freqT = ZIN0(2);
    
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    krC = std::round(freqC / fond0);
    krT = std::round(freqT / fond0);
    
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        //inverse et transpose
        m_k = (float)(std::abs( krC * krT / (i + 1.f) ));
        
        //ajustements
        m_k = std::round( numbins / m_k );
        //m_k = std::max( m_k , 1.f ); //sécurité
        m_k = std::round( numbins / m_k);
        //m_k = std::max( m_k , 1.f ); //sécurité
        
        //réassigne
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Inverse_0_Ctor(PV_Inverse_0 *unit)
{
    SETCALC(PV_Inverse_0_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Inverse_0_Dtor(PV_Inverse_0 *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}


///////////////////////////////////////////////////////////////////////////////////////////////
//*************
//PV_Inverse_1 (original safe)
//*************

void PV_Inverse_1_next(PV_Inverse_1 *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr, freqC, freqT,fond0;
    float krC,krT;
    int size;
    
    freqC = ZIN0(1);
    freqT = ZIN0(2);
    
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    krC = std::round( freqC / fond0 ) ;
    krT = std::round( freqT / fond0 ) ;
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        //inverse et transpose
        m_k = (float)( std::abs( krC * krT / (i + 1.f) ) );
        
        //ajustements
        m_k = std::round( numbins / m_k );
        m_k = std::max( m_k , 1.f );
        m_k = std::round( numbins / m_k);
        m_k = std::max( m_k , 1.f );
        
        //réassigne
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Inverse_1_Ctor(PV_Inverse_1 *unit)
{
    SETCALC(PV_Inverse_1_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Inverse_1_Dtor(PV_Inverse_1 *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//*************
//PV_Inverse_2 (variant safe)
//*************

void PV_Inverse_2_next(PV_Inverse_2 *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr, freqC, freqT,fond0;
    float krC,krT;
    int size;
    
    freqC = ZIN0(1);
    freqT = ZIN0(2);
    
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    krC = freqC / fond0 ;
    krT = freqT / fond0 ;
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        //inverse et transpose
        m_k = (float)( std::abs( krC * krT / (i + 1.f) ) );
        
        //ajustements
        m_k = std::round( numbins / m_k );
        m_k = std::max( m_k , 1.f );
        m_k = std::round( numbins / m_k);
        m_k = std::max( m_k , 1.f );
        
        //réassigne
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Inverse_2_Ctor(PV_Inverse_2 *unit)
{
    SETCALC(PV_Inverse_2_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Inverse_2_Dtor(PV_Inverse_2 *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}


///////////////////////////////////////////////////////////////////////////////////////////////
//*************
//PV_Inverse_Exp (variant safe)
//*************

int moyG (float x){
    return std::sqrt(x * (x + 1.f ) );
};

void PV_Inverse_Exp_next(PV_Inverse_Exp *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr, freqC, freqT,fond0;
    float krC,krT;
    int size;
    
    freqC = ZIN0(1);
    freqT = ZIN0(2);
    
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    krC = std::floor( freqC / fond0 ) ;
    krT = std::floor( freqT / fond0 );
    
    krC = moyG(krC);
    krT = moyG(krT);
    
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        //inverse et transpose
        m_k = (float) (i + 1.f);
        m_k =  moyG(m_k);
        m_k = (float)( std::abs( krC * krT / m_k ) );
        
        //ajustements
        m_k = std::floor( numbins / m_k );
        m_k = moyG(m_k);
        m_k = std::max( m_k , 1.f );
        m_k = std::floor( numbins / m_k);
        m_k = std::max( m_k , 1.f );
        
        //réassigne
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Inverse_Exp_Ctor(PV_Inverse_Exp *unit)
{
    SETCALC(PV_Inverse_Exp_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Inverse_Exp_Dtor(PV_Inverse_Exp *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}

/////////////////////////////////////////////////////////////////////////////////////////////

//******
//PV_Inv
//******

void PV_Inv_next(PV_Inv *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr,fond0;
    int size;
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    float freqC,krC;
    freqC = ZIN0(1);
    
    krC = std::floor(freqC / fond0) + 0.5f;
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        m_k = (float)(krC * krC / (i + 1.5f));

        m_k = std::floor(m_k);
        m_k = std::max( m_k , 1.f );
        
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Inv_Ctor(PV_Inv *unit)
{
    SETCALC(PV_Inv_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Inv_Dtor(PV_Inv *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}


///////////////////////////////////////////////////////////////////////////////////////////////

//*********
//PV_Transp
//*********


void PV_Transp_next(PV_Transp *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float kTransp = ZIN0(1);
    kTransp = std::abs(kTransp);
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        m_k = (float)(kTransp * (i + 1.5f));
        
        m_k = std::floor(m_k);
        m_k = std::max( m_k , 1.f );
     
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_Transp_Ctor(PV_Transp *unit)
{
    SETCALC(PV_Transp_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_Transp_Dtor(PV_Transp *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}



///////////////////////////////////////////////////////////////////////////////////////////////

//*************
//PV_InvTransp
//*************

void PV_InvTransp_next(PV_InvTransp *unit, int inNumSamples)
{
    
    PV_GET_BUF
    MAKE_TEMP_BUF
    
    SCPolarBuf *p = ToPolarApx(buf);
    SCPolarBuf *q = (SCPolarBuf*)unit->m_tempbuf;
    
    // initialize output buf to zeroes
    for (int i=0; i<numbins; ++i) {
        q->bin[i].mag = 0.f;
        q->bin[i].phase = p->bin[i].phase;
    }
    
    float sr,fond0;
    int size;
    
    sr = SAMPLERATE * BUFRATE;
    size = buf->samples;
    fond0 = (float)(sr / size);
    
    
    float freqC, freqT;
    float krC,krT;
    
    freqC = ZIN0(1);
    freqT = ZIN0(2);
    
    krC = std::floor(freqC / fond0) + 0.5 ;
    krT = std::floor(freqT / fond0) + 0.5;
    
    
    q->dc = p->dc;
    q->nyq = p->nyq;
    
    
    float m_k;
    for (int i=0; i < numbins; ++i) {
        
        m_k = (float)std::abs( krC * krT / (i + 1.5f) );        
        m_k = std::floor(m_k);
        m_k = std::max( m_k , 1.f );
        
        int32 iv = (int32)(m_k - 1);
        q->bin[i].mag = p->bin[iv].mag;
        q->bin[i].phase = p->bin[iv].phase;
    };
    
    memcpy(p->bin, q->bin, numbins * sizeof(SCComplex));
    
}

void PV_InvTransp_Ctor(PV_InvTransp *unit)
{
    SETCALC(PV_InvTransp_next);
    ZOUT0(0) = ZIN0(0);
    unit->m_tempbuf = 0;
}

void PV_InvTransp_Dtor(PV_InvTransp *unit)
{
    RTFree(unit->mWorld, unit->m_tempbuf);
}


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void init_SCComplex(InterfaceTable *inTable);

//#define DefinePVUnit(name) \
//    (*ft->fDefineUnit)(#name, sizeof(PV_Unit), (UnitCtorFunc)&name##_Ctor, 0, 0);

#define DefineDtorUnit(name) \
(*ft->fDefineUnit)(#name, sizeof(name), (UnitCtorFunc)&name##_Ctor, \
(UnitDtorFunc)&name##_Dtor, 0);


PluginLoad(PV_Inverse)
{
    ft =inTable;
    init_SCComplex(inTable);
    
    DefineDtorUnit(PV_DistrN);
    
    DefineDtorUnit(PV_Inverse_0);
    DefineDtorUnit(PV_Inverse_1);
    DefineDtorUnit(PV_Inverse_2);
    
    DefineDtorUnit(PV_Inverse_Exp);
    
    DefineDtorUnit(PV_Inv);
    DefineDtorUnit(PV_Transp)
    DefineDtorUnit(PV_InvTransp)
}


