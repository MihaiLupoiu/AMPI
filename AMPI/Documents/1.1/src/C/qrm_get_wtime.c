/* ##############################################################################################
**
** Copyright 2012 CNRS, INPT
**  
** This file is part of qr_mumps.
**  
** qr_mumps is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as 
** published by the Free Software Foundation, either version 3 of 
** the License, or (at your option) any later version.
**  
** qr_mumps is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**  
** You can find a copy of the GNU Lesser General Public License
** in the qr_mumps/doc directory.
**
** ##############################################################################################*/


/*##############################################################################################*/
/** @file qrm_get_wtime.c
 * This file contains timing routines
 *
 * $Date: 2015-10-27 22:41:03 +0100 (mar., 27 oct. 2015) $
 * $Author: abuttari $
 * $Version: 1.1$
 * $Revision: 1980 $
 *
 **/
/*##############################################################################################*/


#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

double qrm_swtime()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return  ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );

}


double qrm_uwtime(){
  struct timeval t;
  struct timezone tzp;

  gettimeofday(&t,&tzp);
  return (double) t.tv_sec*1000000+t.tv_usec;

}



void qrm_msleep(int n){

    struct timespec req={0},rem={0};
    time_t sec=(int)(n/1000);
    n=n-(sec*1000);
    req.tv_sec=sec;
    req.tv_nsec=n*1000000L;
    nanosleep(&req,&rem);
    return;
}
