// MLCG.h                -*- C++ -*-
/* 
Copyright (C) 1988 Free Software Foundation
    written by Dirk Grunwald (grunwald@cs.uiuc.edu)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#ifndef _MLCG_hhh
#define _MLCG_hhh 1 

#include <math.h>
#include "RNG.h"

//
//	Multiplicative Linear Conguential Generator
//

class MLCG : public RNG {
    int initialSeedOne;
    int initialSeedTwo;
    int seedOne;
    int seedTwo;

protected:

public:
    MLCG(int seed1 = 0, int seed2 = 1);
    //
    // Return a long-words word of random bits
    //
    virtual unsigned int asLong();
    virtual void reset();
    int seed1();
    void seed1(int);
    int seed2();
    void seed2(int);
    void reseed(int, int);
};

inline int
MLCG::seed1()
{
    return(seedOne);
}

inline void
MLCG::seed1(int s)
{
    initialSeedOne = s;
    reset();
}

inline int
MLCG::seed2()
{
    return(seedTwo);
}

inline void
MLCG::seed2(int s)
{
    initialSeedTwo = s;
    reset();
}

inline void
MLCG::reseed(int s1, int s2)
{
    initialSeedOne = s1;
    initialSeedTwo = s2;
    reset();
}

#endif
