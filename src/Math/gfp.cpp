/*
Copyright (c) 2017, The University of Bristol, Senate House, Tyndall Avenue, Bristol, BS8 1TH, United Kingdom.
Copyright (c) 2020, COSIC-KU Leuven, Kasteelpark Arenberg 10, bus 2452, B-3001 Leuven-Heverlee, Belgium.

All rights reserved
*/

#include "Math/gfp.h"

#include "Exceptions/Exceptions.h"

template<>
thread_local Zp_Data gfp::ZpD;

template<>
void gfp::AND(const gfp &x, const gfp &y)
{
  bigint bi1, bi2;
  to_bigint(bi1, x);
  to_bigint(bi2, y);
  mpz_and(bi1.get_mpz_t(), bi1.get_mpz_t(), bi2.get_mpz_t());
  to_gfp(*this, bi1);
}

template<>
void gfp::OR(const gfp &x, const gfp &y)
{
  bigint bi1, bi2;
  to_bigint(bi1, x);
  to_bigint(bi2, y);
  mpz_ior(bi1.get_mpz_t(), bi1.get_mpz_t(), bi2.get_mpz_t());
  to_gfp(*this, bi1);
}

template<>
void gfp::XOR(const gfp &x, const gfp &y)
{
  bigint bi1, bi2;
  to_bigint(bi1, x);
  to_bigint(bi2, y);
  mpz_xor(bi1.get_mpz_t(), bi1.get_mpz_t(), bi2.get_mpz_t());
  to_gfp(*this, bi1);
}

template<>
void gfp::AND(const gfp &x, const bigint &y)
{
  bigint bi;
  to_bigint(bi, x);
  mpz_and(bi.get_mpz_t(), bi.get_mpz_t(), y.get_mpz_t());
  to_gfp(*this, bi);
}

template<>
void gfp::OR(const gfp &x, const bigint &y)
{
  bigint bi;
  to_bigint(bi, x);
  mpz_ior(bi.get_mpz_t(), bi.get_mpz_t(), y.get_mpz_t());
  to_gfp(*this, bi);
}

template<>
void gfp::XOR(const gfp &x, const bigint &y)
{
  bigint bi;
  to_bigint(bi, x);
  mpz_xor(bi.get_mpz_t(), bi.get_mpz_t(), y.get_mpz_t());
  to_gfp(*this, bi);
}

template<>
void gfp::SHL(const gfp &x, int n)
{
  if (n < 0)
    {
      throw arithmetic_bug();
    }
  else if (n == 0)
    {
      a= x.a;
    }
  else if (!x.is_zero())
    {
      bigint bi;
      to_bigint(bi, x, false);
      mpn_lshift(bi.get_mpz_t()->_mp_d, bi.get_mpz_t()->_mp_d,
                 bi.get_mpz_t()->_mp_size, n);
      to_gfp(*this, bi);
    }
  else
    {
      assign_zero();
    }
}

template<>
void gfp::SHR(const gfp &x, int n)
{
  if (n < 0)
    {
      throw arithmetic_bug();
    }
  else if (n == 0)
    {
      a= x.a;
    }
  else if (!x.is_zero())
    {
      bigint bi;
      to_bigint(bi, x);
      mpn_rshift(bi.get_mpz_t()->_mp_d, bi.get_mpz_t()->_mp_d,
                 bi.get_mpz_t()->_mp_size, n);
      to_gfp(*this, bi);
    }
  else
    {
      assign_zero();
    }
}

template<>
void gfp::SHL(const gfp &x, const bigint &n)
{
  SHL(x, mpz_get_si(n.get_mpz_t()));
}

template<>
void gfp::SHR(const gfp &x, const bigint &n)
{
  SHR(x, mpz_get_si(n.get_mpz_t()));
}

template<>
gfp gfp::sqrRoot()
{
  // Temp move to bigint so as to call sqrRootMod
  bigint ti;
  to_bigint(ti, *this);
  ti= sqrRootMod(ti, ZpD.pr);
  if (!isOdd(ti))
    ti= ZpD.pr - ti;
  gfp temp;
  to_gfp(temp, ti);
  return temp;
}
