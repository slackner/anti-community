/* Program to do group selection according to the method of Riolo,
 * Cantwell, Reinert, and Newman
 *
 * Written by Mark Newman, January 2017
 */

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "readgml.h"

/* Constants */

#define K 25               // Maximum number of groups
#define MCSWEEPS 10000     // Number of Monte Carlo sweeps to perform
#define SAMPLE 1           // Rate at which to print out results, in sweeps

/* Globals */

NETWORK G;            // The network
int twom;             // Twice the number of edges
double p;             // Average edge probability

int k;                // Current value of k
int *g;               // Group assignments
int *n;               // Group sizes
int **m;              // Edge counts
int *kappa;           // Sums of degrees
int **in;             // Lists of the nodes in each group

double *lnfact;       // Look-up table of log-factorials
gsl_rng *rng;         // Random number generator


// Make a lookup table of log-factorials

void maketable()
{
  int t;
  int length;

  length = twom + G.nvertices + 1;
  lnfact = malloc(length*sizeof(double));
  for (t=0; t<length; t++) lnfact[t] = gsl_sf_lnfact(t);
}


// Initial group assignment

void initgroups()
{
  int i,j,u,v;
  int r;
  int *perm;

  g = malloc(G.nvertices*sizeof(int));

  // Place one randomly chosen node in each group, then put the remaining
  // nodes in random groups.  Works by generating a random permutation of
  // the integers 1...n then placing nodes in groups in that order.

  perm = malloc(G.nvertices*sizeof(int));
  for (i=0; i<G.nvertices; i++) {
    j = gsl_rng_uniform_int(rng,i+1);
    perm[i] = perm[j];
    perm[j] = i;
  }
  for (i=0; i<K; i++) g[perm[i]] = i;
  for (; i<G.nvertices; i++) g[perm[i]] = gsl_rng_uniform_int(rng,K);
  free(perm);

  // Calculate the n's, kappa's, and the lists of group members

  n = calloc(K,sizeof(int));
  kappa = calloc(K,sizeof(int));
  in = malloc(K*sizeof(int*));
  for (r=0; r<K; r++) in[r] = malloc(G.nvertices*sizeof(int));
  for (u=0; u<G.nvertices; u++) {
    r = g[u];
    in[r][n[r]] = u;
    n[r] += 1;
    kappa[r] += G.vertex[u].degree;
  }

  // Calcalate the m's

  m = malloc(K*sizeof(int*));
  for (r=0; r<K; r++) m[r] = calloc(K,sizeof(int));
  for (u=0; u<G.nvertices; u++) {
    for (i=0; i<G.vertex[u].degree; i++) {
      v = G.vertex[u].edge[i].target;
      if ((g[u]!=g[v])||(u<v)) m[g[u]][g[v]]++;
    }
  }

  // Initialize k

  k = K;
}


// Function to calculate the change in log(P(A|g,k)) for a type 1 move

double deltalogp1(int u, int d[], int r, int s)
{
  int t;
  int nrp,nsp,kapparp,kappasp;
  double res;

  if (n[r]==1) {

    // Calculate the new n_s and kappa_s

    nsp = n[s] + 1;
    kappasp = kappa[s] + G.vertex[u].degree;

    // Degree-correction terms

    res = lnfact[G.vertex[u].degree]
           + kappasp*log(nsp) - kappa[s]*log(n[s])
           + lnfact[nsp-1] - lnfact[n[s]-1]
           - lnfact[kappasp+nsp-1] + lnfact[kappa[s]+n[s]-1];

    // Terms in m

    for (t=0; t<k; t++) {
      if (t==r) {
	res += log(0.5*p+1);
      } else if (t==s) {
        res += lnfact[m[s][s]+d[s]] - lnfact[m[s][s]]
                 - (m[s][s]+d[s]+1)*log(0.5*p*nsp*nsp+1)
                 + (m[s][s]+1)*log(0.5*p*n[s]*n[s]+1);
      } else {
	res += lnfact[m[s][t]+d[t]] - lnfact[m[s][t]]
                 - (m[s][t]+d[t]+1)*log(p*nsp*n[t]+1)
                 + (m[s][t]+1)*log(p*n[s]*n[t]+1)
	         - lnfact[d[t]] + (d[t]+1)*log(p*n[t]+1);
      }
    }
    res -= lnfact[d[s]] - (d[s]+1)*log(p*n[s]+1);

  } else {

    // Calculate the new n's and kappa's

    nrp = n[r] - 1;
    nsp = n[s] + 1;
    kapparp = kappa[r] - G.vertex[u].degree;
    kappasp = kappa[s] + G.vertex[u].degree;

    // Degree-correction terms

    res = kapparp*log(nrp) - kappa[r]*log(n[r])
           + lnfact[nrp-1] - lnfact[n[r]-1]
           - lnfact[kapparp+nrp-1] + lnfact[kappa[r]+n[r]-1]
           + kappasp*log(nsp) - kappa[s]*log(n[s])
           + lnfact[nsp-1] - lnfact[n[s]-1]
           - lnfact[kappasp+nsp-1] + lnfact[kappa[s]+n[s]-1];

    // Terms in m

    for (t=0; t<k; t++) {
      if (t==r) {
	res += lnfact[m[r][r]-d[r]] - lnfact[m[r][r]]
                 - (m[r][r]-d[r]+1)*log(0.5*p*nrp*nrp+1)
                 + (m[r][r]+1)*log(0.5*p*n[r]*n[r]+1);
      } else if (t==s) {
	res += lnfact[m[s][s]+d[s]] - lnfact[m[s][s]]
                 - (m[s][s]+d[s]+1)*log(0.5*p*nsp*nsp+1)
                 + (m[s][s]+1)*log(0.5*p*n[s]*n[s]+1);
      } else {
	res += lnfact[m[r][t]-d[t]] - lnfact[m[r][t]]
                 - (m[r][t]-d[t]+1)*log(p*nrp*n[t]+1)
                 + (m[r][t]+1)*log(p*n[r]*n[t]+1)
	         + lnfact[m[s][t]+d[t]] - lnfact[m[s][t]]
                 - (m[s][t]+d[t]+1)*log(p*nsp*n[t]+1)
                 + (m[s][t]+1)*log(p*n[s]*n[t]+1);
      }
    }
    res += lnfact[m[r][s]+d[r]-d[s]] - lnfact[m[r][s]]
             - (m[r][s]+d[r]-d[s]+1)*log(p*nrp*nsp+1)
             + (m[r][s]+1)*log(p*n[r]*n[s]+1);

  }

  return res;
}


// Function to calculate the change in log(P(A|g,k)) for a type 2 move

double deltalogp2(int u, int d[], int r)
{
  int t;
  int nrp,kapparp;
  double res;

  // Calculate the new n_r and kappa_r

  nrp = n[r] - 1;
  kapparp = kappa[r] - G.vertex[u].degree;

  // Degree-correction terms

  res = kapparp*log(nrp) - kappa[r]*log(n[r])
          + lnfact[nrp-1] - lnfact[n[r]-1]
          - lnfact[kapparp+nrp-1] + lnfact[kappa[r]+n[r]-1]
          - lnfact[G.vertex[u].degree];

  // Terms in m

  for (t=0; t<k; t++) {
    if (t==r) {
      res += lnfact[m[r][r]-d[r]] - lnfact[m[r][r]]
               - (m[r][r]-d[r]+1)*log(0.5*p*nrp*nrp+1)
               + (m[r][r]+1)*log(0.5*p*n[r]*n[r]+1);
    } else {
      res += lnfact[m[r][t]-d[t]] - lnfact[m[r][t]]
               - (m[r][t]-d[t]+1)*log(p*nrp*n[t]+1)
               + (m[r][t]+1)*log(p*n[r]*n[t]+1)
	       + lnfact[d[t]] - (d[t]+1)*log(p*n[t]+1);
    }
  }
  res += lnfact[d[r]] - (d[r]+1)*log(p*nrp+1) - log(0.5*p+1);

  return res;
}


// Function to calculate the complete log probability, up to constants

double energy()
{
  int r,s;
  double res;

  res = -k*log(G.nvertices-1);
  for (r=0; r<k; r++) {
    res += lnfact[n[r]];
    res += kappa[r]*log(n[r]) + lnfact[n[r]-1] - lnfact[kappa[r]+n[r]-1];
    res += lnfact[m[r][r]] - (m[r][r]+1)*log(0.5*p*n[r]*n[r]+1);
    for (s=r+1; s<k; s++) {
      res += lnfact[m[r][s]] - (m[r][s]+1)*log(p*n[r]*n[s]+1);
    }
  }

  return res;
}


// Function to update n, kappa, and m for a proposed move

void nmupdate(int u, int d[], int r, int s)
{
  int t;

  n[r]--;
  n[s]++;
  kappa[r] -= G.vertex[u].degree;
  kappa[s] += G.vertex[u].degree;
  for (t=0; t<k; t++) {
    if (t==r) {
      m[r][r] -= d[r];
      m[s][r] += d[r];
      m[r][s] += d[r];
    } else if (t==s) {
      m[r][s] -= d[s];
      m[s][r] -= d[s];
      m[s][s] += d[s];
    } else {
      m[r][t] -= d[t];
      m[t][r] -= d[t];
      m[s][t] += d[t];
      m[t][s] += d[t];
    }
  }
}


// Function that does one Monte Carlo sweep (i.e., n individual moves)

double sweep()
{
  int i,j,u,v;
  int r,s,t;
  int qr;
  int accept=0;
  int d[K];

  for (i=0; i<G.nvertices; i++) {

    // Choose the type of move we are going to do

    if (gsl_rng_uniform_int(rng,G.nvertices-1)>0) {

      // Type 1: move a node

      if (k==1) continue;         // No moves are possible when k=1

      // Choose two distinct groups at random

      r = gsl_rng_uniform_int(rng,k);
      s = gsl_rng_uniform_int(rng,k-1);
      if (r==s) s = k - 1;

      // Choose a node u at random from group r

      qr = gsl_rng_uniform_int(rng,n[r]);
      u = in[r][qr];

      // Count the number of edges this node has to each group

      for (t=0; t<k; t++) d[t] = 0;
      for (j=0; j<G.vertex[u].degree; j++) {
	v = G.vertex[u].edge[j].target;
	d[g[v]]++;
      }

      // Decide whether to accept the move

      if (gsl_rng_uniform(rng)<exp(deltalogp1(u,d,r,s))) {
	nmupdate(u,d,r,s);
	g[u] = s;
	in[r][qr] = in[r][n[r]];
	in[s][n[s]-1] = u;
	accept++;

	// If we have removed the last node from group r,
        // decrease k and relabel

	if (n[r]==0) {

	  k -= 1;
	  if (r!=k) {

	    // Update the group labels and the lists

	    for (qr=0; qr<n[k]; qr++) {
	      v = in[r][qr] = in[k][qr];
	      g[v] = r;
	    }

	    // Update n_r and kappa_r

	    n[r] = n[k];
	    n[k] = 0;
	    kappa[r] = kappa[k];
	    kappa[k] = 0;

	    // Update m_rs
	  
	    for (t=0; t<k; t++) {
	      if (t==r) {
		m[r][r] = m[k][k];
		m[k][k] = 0;
	      } else {
		m[r][t] = m[k][t];
		m[t][r] = m[t][k];
		m[k][t] = m[t][k] = 0;
	      }
	    }
	  }
	}
      }

    } else {

      // Type 2: move a node to a newly created group

      if (k==K) continue;       // No room to increase k, so do nothing

      // Choose source group at random

      r = gsl_rng_uniform_int(rng,k);
      if (n[r]==1) continue;    // Moving last node in group does nothing

      // Choose a node at random in group r

      qr = gsl_rng_uniform_int(rng,n[r]);
      u = in[r][qr];

      // Count the number of edges this node has to each group

      for (t=0; t<k+1; t++) d[t] = 0;
      for (j=0; j<G.vertex[u].degree; j++) {
	v = G.vertex[u].edge[j].target;
	d[g[v]]++;
      }
      
      // Decide whether to accept the move

      if (gsl_rng_uniform(rng)<exp(deltalogp2(u,d,r))) {
	k += 1;
        nmupdate(u,d,r,k-1);
	g[u] = k - 1;
	in[r][qr] = in[r][n[r]];
	in[k-1][0] = u;
	accept++;
      }
    }
  }

  return (double)accept/G.nvertices;
}


main(int argc, char *argv[])
{
  int u,r,s;
  double pr,entropy;
  FILE *f;

  // Initialize the random number generator from the system clock

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,time(NULL));

  // Read the network from stdin

  read_network(&G,stdin);
  for (u=twom=0; u<G.nvertices; u++) twom += G.vertex[u].degree;
  p = (double)twom/(G.nvertices*G.nvertices);
  fprintf(stderr,"Read network with %i nodes and %i edges and p = %g\n",
	  G.nvertices,twom/2,p);

  // Make lookup table of log-factorials

  maketable();

  // Initialize the group assignment

  initgroups();

  // Perform the Monte Carlo

  for (s=0; s<MCSWEEPS; s++) {
    sweep();
    if (s%SAMPLE==0) {

      // Calculate entropy of group assignment

      entropy = 0.0;
      for (r=0; r<k; r++) {
	pr = (double)n[r]/G.nvertices;
	entropy -= pr*log(pr);
      }

      // Print out results

      printf("%i %i %g %g\n",s,k,exp(entropy),energy());
    }
  }
}
