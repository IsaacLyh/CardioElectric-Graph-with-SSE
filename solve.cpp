/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "Plotting.h"
#include  <emmintrin.h>

#ifdef SSE_VEC
#include <emmintrin.h>
#endif 
using namespace std;

void repNorms(double l2norm, double mx, double dt, int m,int n, int niter, int stats_freq);
void stats(double *E, int m, int n, double *_mx, double *sumSq);

#ifdef SSE_VEC
// If you intend to vectorize using SSE instructions, you must
// disable the compiler's auto-vectorizer
__attribute__((optimize("no-tree-vectorize")))
#endif 


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
double L2Norm(double sumSq){
    double l2norm = sumSq /  (double) ((cb.m)*(cb.n));
    l2norm = sqrt(l2norm);
    return l2norm;
}

void solve(double **_E, double **_E_prev, double *R, double alpha, double dt, Plotter *plotter, double &L2, double &Linf){

 // Simulated time is different from the integer timestep number
 double t = 0.0;

 double *E = *_E, *E_prev = *_E_prev;
 double *R_tmp = R;
 double *E_tmp = *_E;
 double *E_prev_tmp = *_E_prev;
 double mx, sumSq;
 int niter;
 int m = cb.m, n=cb.n;
 int innerBlockRowStartIndex = (n+2)+1;
 int innerBlockRowEndIndex = (((m+2)*(n+2) - 1) - (n)) - (n+2);


 // We continue to sweep over the mesh until the simulation has reached
 // the desired number of iterations
  for (niter = 0; niter < cb.niters; niter++){
      if  (cb.debug && (niter==0)){
	  stats(E_prev,m,n,&mx,&sumSq);
          double l2norm = L2Norm(sumSq);
	  repNorms(l2norm,mx,dt,m,n,-1, cb.stats_freq);
	  if (cb.plot_freq)
	      plotter->updatePlot(E,  -1, m+1, n+1);
      }

   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    * These are physical boundary conditions, and are not to be confused
    * with ghost cells that we would use in an MPI implementation
    *
    * The reason why we copy boundary conditions is to avoid
    * computing single sided differences at the boundaries
    * which increase the running time of solve()
    *
    */
    
    // 4 FOR LOOPS set up the padding needed for the boundary conditions
    int i,j;

    // Fills in the TOP Ghost Cells
    for (i = 0; i < (n+2); i++) {
        E_prev[i] = E_prev[i + (n+2)*2];
    }

    // Fills in the RIGHT Ghost Cells
    for (i = (n+1); i < (m+2)*(n+2); i+=(n+2)) {
        E_prev[i] = E_prev[i-2];
    }

    // Fills in the LEFT Ghost Cells
    for (i = 0; i < (m+2)*(n+2); i+=(n+2)) {
        E_prev[i] = E_prev[i+2];
    }	

    // Fills in the BOTTOM Ghost Cells
    for (i = ((m+2)*(n+2)-(n+2)); i < (m+2)*(n+2); i++) {
        E_prev[i] = E_prev[i - (n+2)*2];
    }

    //initialize the SSE BASESTRUCTURE
    __m128d vec1, vec2, vec3, vec4;
    __m128d alpha_128 = _mm_set_pd1(alpha);
    __m128d four_128 = _mm_set_pd1((double)4);
    
    // Solve for the excitation, a PDE
    for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
        E_tmp = E + j;
        E_prev_tmp = E_prev + j;
        //E_tmp[0] = E_prev_tmp[0]+alpha*(E_prev_tmp[1]+E_prev_tmp[-1]-4*E_prev_tmp[0]+E_prev_tmp[0+(n+2)]+E_prev_tmp[0-(n+2)]);
        //E_tmp[1] = E_prev_tmp[1]+alpha*(E_prev_tmp[2]+E_prev_tmp[1-1]-4*E_prev_tmp[1]+E_prev_tmp[1+(n+2)]+E_prev_tmp[1-(n+2)]);
        for(i = 0; i < n; i+=2){
          //load first three argu;
            vec1 = _mm_loadu_pd(&E_prev_tmp[i]); 
            vec2 = _mm_load_pd(&E_prev_tmp[i+1]);
            vec3 = _mm_load_pd(&E_prev_tmp[i-1]);
          //apply multiplication
            vec2 = _mm_mul_pd(vec2,alpha_128);
            vec3 = _mm_mul_pd(vec3,alpha_128);

            vec4 = _mm_add_pd(vec1,vec2);
            vec4 = _mm_add_pd(vec4,vec3);
          //compute forthval
            vec1 = _mm_mul_pd(vec1,alpha_128);
            vec1 = _mm_mul_pd(vec1,four_128);
            vec4 = _mm_sub_pd(vec4,vec1);
          //compute last 2 arguments;
            vec1 = _mm_loadu_pd(&E_prev_tmp[i+(n+2)]);
            vec2 = _mm_loadu_pd(&E_prev_tmp[i-(n+2)]);
            vec1 = _mm_mul_pd(vec1,alpha_128);
            vec2 = _mm_mul_pd(vec2,alpha_128);
            vec4 = _mm_add_pd(vec4,vec1);
            vec4 = _mm_add_pd(vec4,vec2);
            _mm_storeu_pd(&E_tmp[i],vec4);
            //E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
        }
        //E_tmp[n-2] = E_prev_tmp[n-2]+alpha*(E_prev_tmp[n-1]+E_prev_tmp[n-3]-4*E_prev_tmp[n-2]+E_prev_tmp[n-2+(n+2)]+E_prev_tmp[n-2-(n+2)]);
        //E_tmp[n-1] = E_prev_tmp[n-1]+alpha*(E_prev_tmp[n]+E_prev_tmp[n-2]-4*E_prev_tmp[n-1]+E_prev_tmp[n-1+(n+2)]+E_prev_tmp[n-1-(n+2)]);
  
    }

    /* 
     * Solve the ODE, advancing excitation and recovery variables
     *     to the next timtestep
     */
    //e_tmp -- dt, kk,a --prefab for E_tmp;
    __m128d dt2,kk2,a2,one2;
    dt2 = _mm_set_pd1(-dt);
    kk2 = _mm_set_pd1(kk);
    a2 = _mm_set_pd1(a);
    one2 = _mm_set_pd1(1);
    //done building prefab for E_tmp;
    //R_tmp, dt,epsilon, M1,M2,kk,b,1
    __m128d dt22,M22,epsilon2,M12,empty2,b2;
    dt22 = _mm_set_pd1(dt);
    b2 = _mm_set_pd1(b);
    empty2 = _mm_set_pd1(0);
    M22  = _mm_set_pd1(M2);
    epsilon2 = _mm_set_pd1(epsilon);
    M12 = _mm_set_pd1(M1);
    //done setting prefab
    for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
        E_tmp = E + j;
        R_tmp = R + j;
        //because the first position of mesh starts on 0, it is not aligned in memeory;
        //thus we compute it first;
        E_tmp[0] += -dt*(kk*E_tmp[0]*(E_tmp[0]-a)*(E_tmp[0]-1)+E_tmp[0]*R_tmp[0]);
        for(i = 1; i < n-1; i+=2) {
            //this one load value into vector for computation;
            vec4 = _mm_load_pd(&E_tmp[i]);
            //apply first part of array;
            vec2 = _mm_sub_pd(vec4,a2);
            vec3 = _mm_sub_pd(vec4,one2);
            //multiple the required argument as in the function;
            vec1 = _mm_mul_pd(kk2,vec4);
            vec1 = _mm_mul_pd(vec1,vec2);
            vec1 = _mm_mul_pd(vec1,vec3);
            //this one further load more data 
            vec3 = _mm_load_pd(&R_tmp[i]);
            vec2 = _mm_mul_pd(vec4,vec3);
            vec1 = _mm_add_pd(vec1,vec2);
            vec1 = _mm_mul_pd(vec1,dt2);
            vec1 = _mm_add_pd(vec4,vec1);
            //store the computed value
            _mm_store_pd(&E_tmp[i],vec1);
            //E_tmp[i] += -dt*(kk*E_tmp[i]*(E_tmp[i]-a)*(E_tmp[i]-1)+E_tmp[i]*R_tmp[i]);
        }
        //these two are used to avoid possible unaligned data
        E_tmp[n-1] += -dt*(kk*E_tmp[n-1]*(E_tmp[n-1]-a)*(E_tmp[n-1]-1)+E_tmp[n-1]*R_tmp[n-1]);
        R_tmp[0] += dt*(epsilon+M1* R_tmp[0]/( E_tmp[0]+M2))*(-R_tmp[0]-kk*E_tmp[0]*(E_tmp[0]-b-1));
        __m128d E_tmp_vec, R_tmp_vec;
        for(i = 1; i < n-1; i+=2) {
            R_tmp_vec = _mm_load_pd(&R_tmp[i]);
            E_tmp_vec = _mm_load_pd(&E_tmp[i]);
            vec2 = _mm_sub_pd(empty2,R_tmp_vec);
            //-R_tmp[i] == vec2
            vec3 = _mm_sub_pd(E_tmp_vec,b2);
            vec3 = _mm_sub_pd(vec3,one2);
            //E_tmp - b - 1 = vec2
            vec4 = _mm_mul_pd(E_tmp_vec,kk2);
            vec4 = _mm_mul_pd(vec4,vec3);
            //kk*E_tmp[i]*(E_tmp[i]-b-1)) = vec4
            vec2 = _mm_sub_pd(vec2,vec4);
            //(-R_tmp[i]-kk*E_tmp[i]*(E_tmp[i]-b-1)) = vec2;
            vec1 = _mm_add_pd(E_tmp_vec,M22);
            //E_tmp[i] + M2 = vec1
            vec3 = _mm_mul_pd(R_tmp_vec,M12);
            vec3 = _mm_div_pd(vec3,vec1);
            vec3 = _mm_add_pd(epsilon2,vec3);
            vec2 = _mm_mul_pd(vec3,vec2);
            vec2 = _mm_mul_pd(dt22,vec2);
            vec2 = _mm_add_pd(R_tmp_vec,vec2);
            _mm_store_pd(&R_tmp[i],vec2);
        }
        //to avoid unligned potions.
        R_tmp[n-1] += dt*(epsilon+M1* R_tmp[n-1]/( E_tmp[n-1]+M2))*(-R_tmp[n-1]-kk*E_tmp[n-1]*(E_tmp[n-1]-b-1));
    }
  
   if (cb.stats_freq){
     if ( !(niter % cb.stats_freq)){
        stats(E,m,n,&mx,&sumSq);
        double l2norm = L2Norm(sumSq);
        repNorms(l2norm,mx,dt,m,n,niter, cb.stats_freq);
    }
   }

   if (cb.plot_freq){
          if (!(niter % cb.plot_freq)){
	    plotter->updatePlot(E,  niter, m, n);
        }
    }

   // Swap current and previous meshes
    double *tmp = E; E = E_prev; E_prev = tmp;

 } //end of 'niter' loop at the beginning

  // return the L2 and infinity norms via in-out parameters
  stats(E_prev,m,n,&Linf,&sumSq);
  L2 = L2Norm(sumSq);

  // Swap pointers so we can re-use the arrays
  *_E = E;
  *_E_prev = E_prev;
}
