// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MTRand
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.1  28 September 2009  wagnerr@umich.edu

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2009, Richard J. Wagner
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote
//      products derived from this software without specific prior written
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
//
//     When you use this, send an email to: m-mat@math.sci.hiroshima-u.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: wagnerr@umich.edu and Cokus@math.washington.edu
// when you write.

#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma gcc diagnostic ignored "-Wdeprecated-register"

#ifndef MTRAND_H
#define MTRAND_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MTRand object)

#include <climits>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>

class MTRand {
  // Data
public:
  typedef unsigned long uint32; // unsigned integer type, at least 32 bits

  enum { N = 624 };      // length of state vector
  enum { SAVE = N + 1 }; // length of array for save()

protected:
  enum { M = 397 }; // period parameter

  uint32 state[N]; // internal state
  uint32 *pNext;   // next value to get from state
  int left;        // number of values left before reload needed

  // Methods
public:
  MTRand(const uint32 oneSeed); // initialize with a simple uint32
  MTRand(uint32 *const bigSeed, uint32 const seedLength = N); // or array
  MTRand(); // auto-initialize with /dev/urandom or time() and clock()
  MTRand(const MTRand &o); // copy

  // Do NOT use for CRYPTOGRAPHY without securely hashing several returned
  // values together, otherwise the generator state can be learned after
  // reading 624 consecutive values.

  // Access to 32-bit random numbers
  uint32 randInt();               // integer in [0,2^32-1]
  uint32 randInt(const uint32 n); // integer in [0,n] for n < 2^32
  // double rand();                        // real number in [0,1] -- disabled
  // by rtoso double rand( const double n );        // real number in [0,n]
  // -- disabled by rtoso
  double randExc();                  // real number in [0,1)
  double randExc(const double n);    // real number in [0,n)
  double randDblExc();               // real number in (0,1)
  double randDblExc(const double n); // real number in (0,n)
  double operator()();               // same as rand53()

  // Access to 53-bit random numbers (capacity of IEEE double precision)
  double rand();   // calls rand53() -- modified by rtoso
  double rand53(); // real number in [0,1)

  // Access to nonuniform random number distributions
  double randNorm(const double mean = 0.0, const double stddev = 1.0);

  // Re-seeding functions with same behavior as initializers
  void seed(const uint32 oneSeed);
  void seed(uint32 *const bigSeed, const uint32 seedLength = N);
  void seed();

  // Saving and loading generator state
  void save(uint32 *saveArray) const; // to array of size SAVE
  void load(uint32 *const loadArray); // from such array
  friend std::ostream &operator<<(std::ostream &os, const MTRand &mtrand);
  friend std::istream &operator>>(std::istream &is, MTRand &mtrand);
  MTRand &operator=(const MTRand &o);

protected:
  void initialize(const uint32 oneSeed);
  void reload();
  uint32 hiBit(const uint32 u) const { return u & 0x80000000UL; }
  uint32 loBit(const uint32 u) const { return u & 0x00000001UL; }
  uint32 loBits(const uint32 u) const { return u & 0x7fffffffUL; }
  uint32 mixBits(const uint32 u, const uint32 v) const {
    return hiBit(u) | loBits(v);
  }
  uint32 magic(const uint32 u) const { return loBit(u) ? 0x9908b0dfUL : 0x0UL; }
  uint32 twist(const uint32 m, const uint32 s0, const uint32 s1) const {
    return m ^ (mixBits(s0, s1) >> 1) ^ magic(s1);
  }
  static uint32 hash(time_t t, clock_t c);
};

#endif
