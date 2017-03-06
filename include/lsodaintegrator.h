/**********************************************************************
QTAIM - Extension for Quantum Theory of Atoms In Molecules Analysis

Copyright (C) 2010 Eric C. Brown

This file is part of the Avogadro molecular editor project.
For more information, see <http://avogadro.openmolecules.net/>
Avogadro is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
Avogadro is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

**********************************************************************/

#ifndef lsodaintegrator_h
#define lsodaintegrator_h


#include <stdio.h>
#include <math.h>

/*

	The LSODA integrator is now implemented as an auxiliary class for iWQModel
	It was declared as friend for any kinds of iWQModels and the main reason 
	for its separate existence as a class is that its numerous static variables 
	can be separated from real iWQModel contents. 
	
	A possible development would be to implement it as a singleton class serving
	multiple models for multiple timesteps. Now each integrator is "one-shot" so
	it is created to make a single timestep and deleted afterwards.
	
*/
class iWQModel;

class iWQLSODAIntegrator
{
	public:
		//an interface to LSODA for iWQModels
		bool solve1Step(iWQModel * m, double tstart, double tend, double eps);
	private:
		// Associated iWQModel instance
		iWQModel * model;
		
		// LSODA integrator
		void f( int neq, double t, double *y, double *ydot );
		void daxpy( int n, double da, double *dx, int incx, double *dy, int incy );
		double ddot( int n, double *dx, int incx, double *dy, int incy );
		void dgefa( double **a, int n, int *ipvt, int *info );
		void dgesl( double **a, int n, int *ipvt, double *b, int job );
		void dscal( int n, double da, double *dx, int incx );
		int idamax( int n, double *dx, int incx );
		void terminate( int *istate );
		void terminate2( double *y, double *t );
		void successreturn( double *y, double *t, int itask, int ihit, double tcrit, int *istate );
		void lsoda( int neq, double *y, double *t, double tout, int itol, double *rtol, double *atol, int itask, int *istate, int iopt, int jt, int iwork1, int iwork2, int iwork5, int iwork6, int iwork7, int iwork8, int iwork9, double rwork1, double rwork5, double rwork6, double rwork7 );
		void stoda( int neq, double *y );
		void ewset( int itol, double *rtol, double *atol, double *ycur );
		void intdy( double t, int k, double *dky, int *iflag );
		void cfode( int meth );
		void scaleh( double  *rh, double *pdh );
		void prja( int neq, double *y );
		double vmnorm( int n, double *v, double *w );
		double fnorm( int n, double **a, double *w );
		void correction( int neq, double *y, int *corflag, double pnorm, double *del, double *delp, double *told, int *ncf, double *rh, int *m );
		void corfailure( double *told, double *rh, int *ncf, int *corflag );
		void solsy( double *y );
		void methodswitch( double dsm, double pnorm, double *pdh, double *rh );
		void endstoda();
		void orderswitch( double *rhup, double dsm, double *pdh, double *rh, int *orderflag );
		void resetcoeff();
		void freevectors();

		/* newly added static variables */
		int ml, mu, imxer;
		int mord[3];
		double sqrteta, *yp1, *yp2;
		double sm1[13];

		/* static variables for lsoda() */
		double ccmax, el0, h, hmin, hmxi, hu, rc, tn;
		int illin, init, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu;
		double tsw, pdnorm;
		int ixpr, jtyp, mused, mxordn, mxords;

		/* no static variable for prja(), solsy() */

		/* static variables for stoda() */
		double conit, crate, el[14], elco[13][14], hold, rmax, tesco[13][4];
		int ialth, ipup, lmax, meo, nslp;
		double pdest, pdlast, ratio, cm1[13], cm2[6];
		int icount, irflag;

		/* static variable for block data */
		int mesflg;

		/* static variables for various vectors and the Jacobian. */
		double **yh, **wm, *ewt, *savf, *acor;
		int *ipvt;
		int m_lenyh;
		int m_nyh;
	};

#endif
