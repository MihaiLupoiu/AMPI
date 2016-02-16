/*
 *  MC_scheduler.c
 *  tiledQR
 *
 *  Created by Mathias Jacquelin on 3/24/11.
 *  Copyright 2011 LIP/ENS-Lyon. All rights reserved.
 *
 */


#include "MC_scheduler.h"
#include "task.h"
#include "common.h"
#include "main_scheduler.h"

void MC_init(){
  
  msched_args_t ps_args = (msched_args_t)MSG_host_get_data(MSG_host_self()); 
  task_t facto,update,elim;

  unsigned int p,q;
  unsigned int i, j, k;

  
  p = ps_args->ui_p;
  q = ps_args->ui_q;
  
  unsigned int r, s, x;
  unsigned int **z;

  
  x = 0; while ( x*(x+1)/2 < p-1 ) x++;
  
  z = (unsigned int**)malloc(p*sizeof(unsigned int *));
  for ( i = 0; i < p; i++) z[i] = (unsigned int*)malloc(q*sizeof(unsigned int));
  
  r = 1;
  for ( s = 0; s < x; s++ ){
    for ( i = r; i < (((r+s+1)<p)?(r+s+1):p); i++){
      z[i][0] = x - s;
    }
    r = r + s + 1;
  }
  
  for ( j = 1; j < q; j++ )
    for ( i = j+1; i < p; i++ )
      z[i][j] = 2 + z[i-1][j-1];
  
  unsigned int t, tmax, inc, first;
  
  
  if ( q < p ) tmax = z[q][q-1];
  if ( q == p ) tmax = z[q-1][q-2];
  
  /* reduce the first column to triangles */
  for (i=0;i<p;i++) {
    
    facto = new_task(F, ps_args,NULL,i, 0,0,0);
    
    for ( j = 1; j < q ; j++ ){
      update = new_task(H,ps_args,facto, i,0, j,0);         
    }
    
    
  }
  
  for ( t = 1; t <= tmax; t++ ){
    for ( j = 0; j < q; j++ ){
      first = 1;
      for ( i = p-1; i > j; i-- ){
        if ( z[i][j] == t ) {
          if ( first ){
            inc = 1; while( z[i-inc][j] == t ) inc++;
            first = 0;
          }
          //                          printf("%%at time %2d -- eliminate A[%2d][%2d] with A[%2d][%2d]\n",t,i,j,i-inc,j);
          /* First we zero out the lower triangle */
          elim = new_task(Z,ps_args,NULL, i, i-inc, j,0);
          
          /* and send the update across the row */
          for (k=j+1;k<q;k++) {
            update = new_task(V,ps_args,elim, i,i-inc, j,k); 
            
          }
          
          
          /* convert right-neighboring tile into triangle */
          if ( j+1 < q ) {
            facto = new_task(F, ps_args,NULL,i, j+1,0,0);
            
            for (k=j+2;k<q;k++) {
              update = new_task(H,ps_args,facto, i,j+1, k,0);         
            }
            
            
          }
        }
      }
    }
  }
  
  free(z);

}
