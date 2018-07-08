/* Program to calculate the number of communities in a network using the
 * method of Newman and Reinert, which calculates a posterior probability
 * by Monte Carlo simulation of the integrated likelihood of a
 * degree-corrected stochastic block model
 *
 * Written by Mark Newman  6 APR 2016
 */

/* Program control */

#define VERBOSE

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "readgml.h"

/* Constants */

#define K 40             // Maximum number of groups
#define MCSWEEPS 10000   // Number of Monte Carlo sweeps
#define SAMPLE 1         // Interval at which to print out results, in sweeps

/* Globals */

NETWORK G;               // Struct storing the network
int twom;                // Twice the number of edges
double p;                // Average edge probability

int k;                   // Current value of k
int *g;                  // Group assignments
int *n;                  // Group sizes
int **m;                 // Edge counts

double *lnfact;          // Look-up table of log-factorials
double E;                // Log probability

gsl_rng *rng;            // Random number generator


// Make a lookup table of log-factorial values

void maketable()
{
  int t;
  int length;

  length = twom + G.nvertices + 1;
  lnfact = malloc(length*sizeof(double));
  for (t=0; t<length; t++) lnfact[t] = gsl_sf_lnfact(t);
}


// Log-probability function

double logp(int *n, int **m)
{
  int r,s;
  int kappa;
  double res=0.0;

  for (r=0; r<k; r++) {
    res += lnfact[n[r]];
    if (n[r]>0) {
      for (s=kappa=0; s<k; s++) kappa += m[r][s];
      res += kappa*log(n[r]) + lnfact[n[r]-1] - lnfact[kappa+n[r]-1];
      res += lnfact[m[r][r]/2] - (m[r][r]/2+1)*log(0.5*p*n[r]*n[r]+1);
      for (s=r+1; s<k; s++) {
    	res += lnfact[m[r][s]] - (m[r][s]+1)*log(p*n[r]*n[s]+1);
      }
    }
  }

  return res;
}


// Initial group assignment

void initgroups()
{
  int i,u,v;
  int r;

  // Make the initial group assignments at random

  g = malloc(G.nvertices*sizeof(int));
  for (u=0; u<G.nvertices; u++) g[u] = gsl_rng_uniform_int(rng,K);

  // Calculate the values of the n's

  n = calloc(K,sizeof(int));
  for (u=0; u<G.nvertices; u++) n[g[u]]++;

  // Calcalate the values of the m's

  m = malloc(K*sizeof(int*));
  for (r=0; r<K; r++) m[r] = calloc(K,sizeof(int));
  for (u=0; u<G.nvertices; u++) {
    for (i=0; i<G.vertex[u].degree; i++) {
      v = G.vertex[u].edge[i].target;
      m[g[u]][g[v]]++;
    }
  }

  // Initialize k and the log-probability

  k = K;
  E = logp(n,m);
}


// Function to update value of k

void changek()
{
  int r,s,u;
  int kp;
  int empty;
  int map[K];
  int sum;

  // With probability 0.5, decrease k, otherwise increase it

  if (gsl_rng_uniform(rng)<0.5) {

    // Count the number of empty groups

    for (r=0,empty=0; r<k; r++) if (n[r]==0) empty++;

    // If there are any empty groups, remove one of them, or otherwise do
    // nothing

    if (empty>0) {

      // If there is more than one empty group, choose at random which one
      // to remove

      do {
        r = gsl_rng_uniform_int(rng,k);
      } while (n[r]>0);

      // Decrease k by 1

      k = k - 1;

      // Update the group labels

      for (u=0; u<G.nvertices; u++) {
	if (g[u]==k) g[u] = r;
      }

      // Update n_r

      n[r] = n[k];

      // Update m_rs

      for (s=0; s<k; s++) {
	if (r==s) {
	  m[r][r] = m[k][k];
	} else {
	  m[r][s] = m[k][s];
	  m[s][r] = m[s][k];
	}
      }
    }

  } else {

    // With probability k/(n+k) increase k by 1, adding an empty group

    if ((G.nvertices+k)*gsl_rng_uniform(rng)<k) {
      if (k<K) {
	n[k] = 0;
	for (r=0; r<=k; r++) m[k][r] = m[r][k] = 0;
	k = k + 1;
      }
    }
  }
}


// Function to update n and m for a proposed move

void nmupdate(int r, int s, int d[])
{
  int t;

  n[r]--;
  n[s]++;
  for (t=0; t<k; t++) {
    m[r][t] -= d[t];
    m[t][r] -= d[t];
    m[s][t] += d[t];
    m[t][s] += d[t];
  }
}


// Function that does one MCMC sweep (i.e., n individual moves) using the
// heatbath algorithm

double sweep()
{
  int i,j,u,v;
  int r,s;
  int temp;
  int accept=0;
  int d[K];
  double x,Z,sum;
  double newE[K];
  double boltzmann[K];

  for (i=0; i<G.nvertices; i++) {

    // Optionally, perform a k-changing move

    if ((G.nvertices+1)*gsl_rng_uniform(rng)<1.0) changek();

    // Choose a random node

    u = gsl_rng_uniform_int(rng,G.nvertices);
    r = g[u];

    // Find the number of edges this node has to each group

    for (s=0; s<k; s++) d[s] = 0;
    for (j=0; j<G.vertex[u].degree; j++) {
      v = G.vertex[u].edge[j].target;
      d[g[v]]++;
    }

    // Calculate the probabilities of moving it to each group in turn

    Z = 0.0;
    for (s=0; s<k; s++) {
      if (s==r) {
	newE[s] = E;
      } else {
	nmupdate(r,s,d);
	newE[s] = logp(n,m);
	nmupdate(s,r,d);
      }
      boltzmann[s] = exp(newE[s]-E);
      Z += boltzmann[s];
    }

    // Choose which move to make based on these probabilities

    x = Z*gsl_rng_uniform(rng);
    for (s=0,sum=0.0; s<k; s++) {
      sum += boltzmann[s];
      if (sum>x) break;
    }

    // Make the move

    if (s!=r) {
      g[u] = s;
      nmupdate(r,s,d);
      E = newE[s];
      accept++;
    }
  }

  return (double)accept/G.nvertices;
}


main(int argc, char *argv[])
{
  int u,r,s;

  // Initialize the random number generator from the system clock

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,time(NULL));

  // Read the network from stdin

#ifdef VERBOSE
  fprintf(stderr,"Reading network...\n");
#endif
  read_network(&G,stdin);
  for (u=twom=0; u<G.nvertices; u++) twom += G.vertex[u].degree;
  p = (double)twom/(G.nvertices*G.nvertices);
#ifdef VERBOSE
  fprintf(stderr,"Read network with %i nodes and %i edges\n",
	  G.nvertices,twom/2);
#endif

  // Make the lookup table

  maketable();

  // Initialize the group assignment

  initgroups();

  // Perform the Monte Carlo

  for (s=0; s<MCSWEEPS; s++) {
    sweep();
    if (s%SAMPLE==0) {
      printf("%i %i %g\n",s,k,E);
#ifdef VERBOSE
      fprintf(stderr,"Sweep %i...\r",s);
#endif
    }
  }
#ifdef VERBOSE
  fprintf(stderr,"\n");
#endif
}
