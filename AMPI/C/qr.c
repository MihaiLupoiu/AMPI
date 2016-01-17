/*
 ============================================================================
 Name        : Assignment AMPI
 Author      : Mihaita Alexandru Lupoiu
 Version     : 0.0.1
 Description : QR Factorization
 ============================================================================
 */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"

#include <string.h>
#include <unistd.h>

/* 
 * Convertir matriz cuadrada entre formato fortran estandar (por columnas,
 * lda=n) y formato de matriz tile con tamaño de bloque bs.
 * A es entrada, B es salida
 * Si tofrom=0, A es la de formato fortran, B es de tipo tile
 * Si tofrom=1, A es la de formato tile, B es de tipo fortran
 */
void convtile(double *A, double *B, int n, int bs, int tofrom)
{
  int i, j, ii, jj, i2, j2, nb=n/bs;

  for (i=0; i<nb; i++) {
    for (j=0; j<nb; j++) {
      for (jj=0; jj<bs; jj++) {
        if (tofrom) 
        	memcpy(B+i*bs+(j*bs+jj)*n, A+(i+j*nb)*bs*bs+jj*bs, bs*sizeof(double));
        else 
        	memcpy(B+(i+j*nb)*bs*bs+jj*bs, A+i*bs+(j*bs+jj)*n, bs*sizeof(double));
      }
    }
  }
}


void QR_LAPACK_Tile(double *A, int lA, int bs, int *info)
{
  	/* The tiled algorithm for QR factorization.
	for k = 1, 2..., min(p, q) do
		DGEQRT(Akk, Vkk, Rkk, Tkk)
		for j = k + 1, k + 2, ..., q do
			DLARFB(Akj , Vkk, Tkk, Rkj )
		end for
		for i = k + 1, k + 1, ..., p do
			DTSQRT(Rkk, Aik, Vik, Tik)
			for j = k + 1, k + 2, ..., q do
				DSSRFB(Rkj , Aij , Vik, Tik)
			end for
		end for
	end for
	*/

  	int i=0, j=0, k=0, nb=lA/bs, lda=bs, m=bs, n=bs, bs2=bs*bs;

  	//DGEQRT(Akk, Vkk, Rkk, Tkk)
  	int ldt=bs;
	double * T = (double *) calloc(ldt * bs, sizeof( double ) );
	double *work = (double *) calloc(n*bs, sizeof( double ) );
	//DLARFB(Akj , Vkk, Tkk, Rkj )
	int K=m, ldc = lda;


	//printf("MATRIZ B\n");

	//printMatrixCustom(A,lA,lA);

	for (k=0; k<nb; k++) { //min(p, q)
		printf("i: %d, j: %d, k: %d, nb: %d, lda: %d, m: %d, n: %d \n",i ,j, k , nb, lda, m, n);

		
		printMatrixCustom(A+(k+k*nb)*bs2,m,n);
		//DGEQRT(Akk, Vkk, Rkk, Tkk)
		//http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational.html#gaddcf152e87deec6123a1899f6f51101e
		dgeqrt_(&m, &n, &bs, A+(k+k*nb)*bs2, &lda, T, &ldt, work, info);
		if (*info != 0) return;
		//printMatrixCustom(A,m,n);

		printMatrixCustom(A+(k+k*nb)*bs2,m,n);
		
		for (j=k+1; j<nb; j++){ //q
			printMatrixCustom(A+(j+k*nb)*bs2,m,n);
			printMatrixCustom(T,bs,bs);
			printf("L, T, m: %d, n: %d, k: %d, nb: %d, V(Akk), ldv(lda): %d, T, ldt: %d, C(Akj), ldc(lda): %d \n", m, n, K, bs, lda, ldt, ldc);

			//DLARFB(Akj , Vkk, Tkk, Rkj )
			//dgemqrt (SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, C, LDC, WORK, INFO)
			//http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational.html#gaf55d7b3137b198647461d429a7e9b2c6

			dgemqrt_("L", "T", &m, &n, &K, &bs, A+(k+k*nb)*bs2, &lda, T, &ldt, A+(j+k*nb)*bs2, &ldc, work, &info);
			if (*info != 0) return;
		}
/*
		for (i=k+1; i<nb; i++){ //p
			DTSQRT(Rkk, Aik, Vik, Tik)
			//http://www.netlib.org/lapack/explore-html/d1/d55/dtpqrt_8f.html#aa02cc2297f978edb5ef2a8fd1dcc9321
			//dtpqrt_(&m, &n, &n, &bs, A, dla, B, ldb, T, ldt, work, info);

			for (j=k+1; j<nb; j++){ //q
				DSSRFB(Rkj , Aij , Vik, Tik)
				//http://www.netlib.org/lapack/explore-html/d7/de6/dtpmqrt_8f.html
				//dtpmqrt_(side, trans, m,n,k,l,nb,v,ldv,t,ldt,A,lad,b,ldb,work,info);
			}
		}
*/
	}

}

void QR_LAPACK(double *A, int lA){
	//http://www.netlib.no/netlib/lapack/double/dgeqrf.f
	/* QR LAPACK */
    int n = lA;
    int m = n;

    int info = 0;
    int k = n;          /* k = min(m,n);       */
    int lda = m;        /* lda = max(m,1);     */
    int lwork = n;      /* lwork = max(n,1);   */
    int max = lwork;    /* max = max(lwork,1); */

    double *work;
    double *tau;
    double *vec;

    work = (double *) calloc(max, sizeof( double ) );
    tau  = (double *) calloc( k, sizeof( double ) );
    vec  = (double *) calloc( m, sizeof( double ) );

    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);

    //printMatrixCustom(A,m,n);
    
    free(work);
    free(tau);
	free(vec);
}

int main(int argc, char *argv[]){

	int major=0, minor=0, patch=0;
    ilaver_(&major,&minor,&patch);
    printf("lapack %d.%d.%d\n",major,minor,patch);


    double *A, *B;
    double *oA;
    int lA = 6, bs = 2, nthreads = 1;
    int info;

    /* Extracción de argumentos */
    if (argc >= 2) {
        if ((lA = atoi(argv[1])) < 0) lA = 6;
    }
    if (argc >= 3) {
        if ((bs = atoi(argv[2])) < 0) bs = 6;
    }
    if (lA<bs || lA%bs) {
        printf("ERROR: el tamaño de la matriz ha de ser divisible por bs\n");
        exit(1);
    }
    if (argc >= 4) {
        if ((nthreads = atoi(argv[3])) < 0) nthreads = 1;
    }
    /* Extracción de argumentos */


    A = (double *) calloc(lA * lA, sizeof(double));
    B = (double *) calloc(lA * lA, sizeof(double));
    oA = (double *) calloc(lA * lA, sizeof(double));

	initializeCust(A,lA);

	copyMatrixCustom(A, oA, lA);
    printMatrixCustom(A,lA,lA);    
    
	/* QR LAPACK */
	
    QR_LAPACK(A,lA);
	printMatrixCustom(A,lA,lA);

	/* QR TILED LAPACK */

	copyMatrixCustom(oA, A, lA);
	//printMatrixCustom(A,lA,lA);
	
	convtile(A, B, lA, bs, 0);
    QR_LAPACK_Tile(B, lA, bs, &info);
    convtile(B, A, lA, bs, 1);

	printMatrixCustom(B,lA,lA);

	/* QR TILED LAPACK */


    free(A);
    free(oA);
    return 0;
}
