#ifndef GSL_H
#define GSL_H

// Routines here are derived from the (G)nu (S)cientific (L)ibrary, version
// 1.5.  Routine-specific copyright notices from the GSL are provided within
// the class definition, though the paragraphs common to each copyright notice
// and to GPL-covered code in general are reproduced below.

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <cmath>

class GSL {
    public:
        GSL();
        ~GSL();

    //////////////////////////////////////  beginning ran3()
    // ran3() random rumber generator, based on Knuth's subtractive algorithm.
    // Modified from Gnu Scientific Library 1.5.  I modified this to work
    // within the GSL class here, with generator state held in a private
    // member.  Separate instantiations of GSL will have separate states.  The
    // copyright notice below is reproduced from the GSL source file.
    /* rng/ran3.c
    *
    * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
    */

    private:
        const long int M_BIG;
        const long int M_SEED;
        bool seed_set;
        struct ran3_state_t {
            unsigned int x;
            unsigned int y;
            unsigned long int buffer[56];
        } state;
        unsigned long int           ran3_get ();
        double                      ran3_get_double ();
        void                        ran3_set (unsigned long int s);
    public:
        // ran3() is my own interface to the private methods that implement the
        // algorithm.  I chose not to inheret the complexity of the GSL
        // interface, but rather come close to the semantics of Numerical
        // Recipes' ran3() as I've implemented them in the past.
        //
        // if seed < 0, set seed using abs(seed)
        // if seed == 0, set seed to default 1
        // if seed > 0 && ! seed_set, then set seed using time(NULL)
        // return random double
        double                      ran3(int seed = 1);

    //
    //////////////////////////////////////  end ran3()

};

//////////////////////////////////////
// constructor with object initialization, and destructor

inline GSL::GSL() :
    // ran3() initializers
    M_BIG(1000000000), M_SEED(161803398), seed_set(false)
{
    /* empty */
}

inline GSL::~GSL() {
    /* empty */
}

//////////////////////////////////////
// beginning of ran3 methods

inline unsigned long int GSL::ran3_get () {
    long int j;
    state.x++;
    if (state.x == 56) 
        state.x = 1;
    state.y++;
    if (state.y == 56) 
        state.y = 1;
    j = state.buffer[state.x] - state.buffer[state.y];
    if (j < 0) 
        j += M_BIG;
    state.buffer[state.x] = j;
    return j;
}

inline double GSL::ran3_get_double () {
    return (ran3_get() / double(M_BIG));
}

inline void GSL::ran3_set (unsigned long int s) {
    int i, i1;
    long int j, k;

    if (s == 0) 
        s = 1;      /* default seed is 1 */
    j = (M_SEED - s) % M_BIG;
    state.buffer[0] = 0;
    state.buffer[55] = j;
    k = 1;
    for (i = 1; i < 55; i++) {
        int n = (21 * i) % 55;
        state.buffer[n] = k;
        k = j - k;
        if (k < 0) 
            k += M_BIG;
        j = state.buffer[n];
    }
    for (i1 = 0; i1 < 4; i1++) {
        for (i = 1; i < 56; i++) {
            long int t = state.buffer[i] - state.buffer[1 + (i + 30) % 55];
            if (t < 0) 
                t += M_BIG;
            state.buffer[i] = t;
        }
    }
    state.x = 0;
    state.y = 31;
    seed_set = true;
    return;
}

inline double GSL::ran3(int seed) {
    if (seed <= 0)
        ran3_set(-seed);
    else if (seed_set == false)
        ran3_set((time( NULL ) % 10000000));
    return(ran3_get_double());
}

// end of ran3 methods
//////////////////////////////////////

#endif // GSL_H


