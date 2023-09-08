#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>


/* 
 * calculates (a * b) % c taking into account that a * b might overflow 
 */

uint64_t mulmod(uint64_t a, uint64_t b, uint64_t mod)
{
  uint64_t x = 0,y = a % mod;
  while (b > 0)
    {
      if (b % 2 == 1)
	{    
	  x = (x + y) % mod;
	}
      y = (y * 2) % mod;
      b /= 2;
    }
  return x % mod;
}

/* 
 * modular exponentiation
 */

uint64_t modulo(uint64_t base, uint64_t exponent, uint64_t mod)
{
  uint64_t x = 1;
  uint64_t y = base;

  while (exponent > 0)
    {
      if (exponent % 2 == 1)
	x = (x * y) % mod;
      y = (y * y) % mod;
      exponent = exponent / 2;
    }

  return x % mod;
}

     

/*
 * Miller-Rabin Primality test, iteration signifies the accuracy
 */

int Miller(uint64_t p,int iteration)
{
  int i;
  uint64_t s;

  if (p < 2)
    {
      return 0;
    }

  if (p != 2 && p % 2==0)
    {
      return 0;
    }

  s = p - 1;
  while (s % 2 == 0)
    {
      s /= 2;
    }

  for (i = 0; i < iteration; i++)
    {
      uint64_t a = rand() % (p - 1) + 1, temp = s;
      uint64_t mod = modulo(a, temp, p);
      while (temp != p - 1 && mod != 1 && mod != p - 1)
	{
	  mod = mulmod(mod, mod, p);
	  temp *= 2;
	}

      if (mod != p - 1 && temp % 2 == 0)
	{
	  return 0;
	}
    }
  return 1;

}


