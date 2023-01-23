/*
Sanjay R Kharche.
*/
/* Write only what you need. */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <unistd.h>

// parameters.
#define DELTAT	0.001
#define Lx 		50.0 // 150.0
#define Ly 		6.4
#define Lxin		10.0 // inlet and outlet of dialysate lengths.
#define Lxout		Lxin // mm. The streamlines do not work unless Lxin == Lxout.
#define dx 		0.2
#define dy 		dx
#define nx 		(int)(floor(Lx/dx))
#define ny 		(int)(floor(Ly/dx))
#define nxdin		(int)floor(Lxin/dx) // this is on top edge to far left.
#define nxdout 	(int)floor(Lxout/dx)

#define aConstant 	0.5*((dx*dx*dy*dy)/(dx*dx+dy*dy)); 

#define Nfibs		75.0      // number of hollow fibres.

#define Qd		16.6 // mm3/s flow of total dialysate at steady state, Janssen et al. sci reps.
#define Qb		8.3 // 1 ml/s or 8.3 mm3/s, Janssen et al. sci reps. flow of total blood at steady state.

// velocities dependent on flow.
// geometry.
#define RM		(Ly/2.0)    	// module inner radius, mm
#define RI      		0.09      		// fiber inner radius, mm.
#define RO      	(RI + 0.045)    // fiber outer radius. Artif Organs, Vol. 33, No. 6, 2009, table 1.
#define L       		Lx 
#define Lwet    	L*1.28          	// Table 1, Chemical Engineerin# Science, Vol. 51, No. 17, pp. 4197-4213, 1996.

#define udl		 (2.0*Qd/(          M_PI*RM*RM)) // mm/s, this comes from Qd, radius.
#define ubr		 -(2.0*Qb/(Nfibs*M_PI*RM*RM)) // mm/s, this comes from Qb, radius.
#define pdR		10.0 // at outlet, dialysate pressure is low.
#define pbL		(67*1333.2) // BIO05001FU.pdf

#define alpha 		0.2
#define betasor	1.2
#define restol		0.05

#define rhod  0.0010295 // 		% g/mm3. Eloot et al. Nephrol Dial Transplant (2007) 22: 2962–2969 doi:10.1093/ndt/gfm356
#define mud  0.0014650 // 		% g/(mm-s). Eloot et al. Nephrol Dial Transplant (2007) 22: 2962–2969 doi:10.1093/ndt/gfm356

#define rhob 0.001050 		// g/mm3. blood density is almost the same as water: https://www.hindawi.com/journals/jvm/2015/152730/
#define mub	0.0031 		// g/(mm-s). see: The International Journal of Artificial Organs / Vol. 32 / no. 6, 2009 / pp. 329-335

// porosities and Darcy permeabilities.
#define epsd 	 (1.0 - Nfibs*RO*RO*Lwet/(RM*RM*L)) // Artif Organs, Vol. 39, No. 6, 2015 has epsd = 0.28, epsb = 0.45.
#define phid 	 (1.0 - epsd)
#define epsb 	 (Nfibs * RI * RI * Lwet / (RM*RM*L))

// dialysate.
#define kdx 		 (RO*RO*L*L/(4.0*phid*Lwet*Lwet*epsd)*(-log(phid) - 3.0/2.0 + 2.0*phid - phid*phid/2.0))
#define kdy 		 (RO*RO/(8.0*phid*epsd)*(-log(phid) + (phid*phid - 1.0)/(phid*phid+1)))
// blood
#define kbx 		 (RI*RI*L*L/(8.0*Lwet*Lwet))
#define kby 		 0.0 // this term must not be written in the formulas for momentum.

// urea.
#define cdlu		0.000175 // g/mm3 is urea concentration in dialysate.
#define cbru		0.001750 // g/mm3. This is urea concentration in blood before dialyser. Fig1 in https://academic.oup.com/ndt/article/16/9/1814/1863029
#define cblu 		cbru/5 // outlet urea, comes from GBH BJ measurements.
#define cddiffu 	0.00182 // mm2/s. this is free diffusion of urea taken from the two porous media paper.
// dialysate diffusions.
#define cddiffux 	cddiffu*(L*L/(Lwet*Lwet))
#define cddiffuy 	cddiffu/(2.0 - epsd) // mm2/s. this is free diffusion of urea taken from the two porous media paper.
// blood diffusions.
#define cbdiffux 	cddiffu*(L*L/(Lwet*Lwet))
#define cbdiffuy 	0.0 // because the hollow fibres are not connected to each other.

// beta is 2.63 mm2/mm3 for now accoring to above paramters.
#define betaa		2.0*Nfibs*RI*Lwet/(RM*RM*L)  // membrane area per unit volume, mm2 /mm3. Eq 7 in porous media paper.

#define Lp		1.15e-7 // mm2-s/g.

// functions.
double residuals(double **now, double **before, int replace){
int i, j, iter;
double res;
res = 0.0;
for(j=0;j<ny;j++)
for(i=0;i<nx;i++){
	res = res + fabs(now[j][i] - before[j][i]);
	if(replace==1) before[j][i] = now[j][i];
	if(!(now[j][i]==now[j][i])){ printf("there is a nan\n"); exit(-2); }
}
return res;
}

void deallocate(double **c, int nny, int nnx){

int i, j;
for(j=0;j<nny;j++) free(c[j]);
free(c);

return;
}


double **doubleAllocate( int nny, int nnx ){
int i, j;
    double **p = ( double ** )calloc( nny , sizeof( double * ) );
    for (j = 0; j < nny; j++ ) 
    p[j] = ( double * )calloc( nnx , sizeof( double ) );
    return p;
}

/* initialize is never called. */
void initialize(double **psyd, double **psyd0){
int i, j;
for(j=0;j<ny;j++){
	psyd0[j][0] 	= udl*Ly*(1.0 - (double)j/(double)(ny-1));
	  psyd[j][0] 	= udl*Ly*(1.0 - (double)j/(double)(ny-1)); // left.
}
for(i=1;i<nx-1;i++){
	psyd0[0][i] 	= udl*Ly;
	  psyd[0][i] 	= udl*Ly;
}
return;
} // end of initialize


void vorticityTransport(double **w, double **w0, double **psy0, double **S1, int whichModel, float DarcyOrNot){

int i, j;
double diff, convec, u, v, darcyterm, S1term;
double rho, mu, kx, ky, dvdx, dudy, eps;
if(whichModel==-1){
rho = rhob; mu = mub;
kx = kbx; ky = 0.0;
eps = epsb;
}
else{
rho = rhod; mu = mud;
kx = kdx; ky = kdy;
eps = epsd;
}

for(j=1;j<ny-1;j++)
for(i=1;i<nx-1;i++){
u 		=    psy0[j+1][i] - psy0[j-1][i]; // u = + d psy / dy.
v 		= -(psy0[j][i+1] - psy0[j][i-1]); // v = -  d psy / dx.
dvdx 	= -(psy0[j][i+1] + psy0[j][i-1] - 2.0*psy0[j][i]);
dudy 	=    psy0[j+1][i] + psy0[j-1][i] - 2.0*psy0[j][i];
diff 		=    w0[j][i+1]    + w0[j][i-1]    + w0[j+1][i] + w0[j-1][i];
convec 	=    u*(w0[j][i+1]-w0[j][i-1]) + v*(w0[j+1][i] - w0[j-1][i]);

w[j][i] 	= 0.25*diff - 0.25*0.25*rho/mu*convec;
}

if(whichModel==-1){
// blood.
	for(j=0;j<ny;j++){
		w[j][0] 		= w[j][1]; // left.
		w[j][nx-1] 	= 0.0; // right.
	}
	for(i=1;i<nx-1;i++){
		w[0][i] 		= 2.0*(psy0[0][i] 		- psy0[1][i])/(dx*dx); // top
		w[ny-1][i] 	= 2.0*(psy0[ny-1][i] 	- psy0[ny-2][i])/(dx*dx);	// bottom
	}
} // end of whichModel = -1, blood.


if(whichModel==1){
// model 1 b.c.
	for(j=0;j<ny;j++){
		w[j][0] 		= 0.0; // left.
		w[j][nx-1] 	= w[j][nx-2]; // right.
	}
	for(i=1;i<nx-1;i++){
		w[0][i] 		= 2.0*(psy0[0][i] 		- psy0[1][i])/(dx*dx); // top
		w[ny-1][i] 	= 2.0*(psy0[ny-1][i] 	- psy0[ny-2][i])/(dx*dx);	// bottom
	}
} // end of whichModel = 1, straight pipe dialysate.

if(whichModel==2){
// model 2 b.c.
	for(j=0;j<ny;j++){
		w[j][0] 		= 2.0*(psy0[j][0] 		- psy0[j][1])/(dx*dx); // left.
		w[j][nx-1] 	= w[j][nx-2]; // right.
	}
	for(i=0;i<nx;i++){
		if(i<nxdin)
			w[0][i] 	= 0.0;
		else
			w[0][i] 	= 2.0*(psy0[0][i] 		- psy0[1][i])/(dx*dx); // top
		
		w[ny-1][i] 	= 2.0*(psy0[ny-1][i] 	- psy0[ny-2][i])/(dx*dx);	// bottom
	}
} // end of whichModel = 2.


if(whichModel==3){
// model 3 b.c.
	for(j=0;j<ny;j++){
		w[j][0] 		= 2.0*(psy0[j][0] 		- psy0[j][1])/(dx*dx); // left.
		w[j][nx-1] 	= w[j][nx-2]; // right.
	}
	for(i=0;i<nx;i++){
		if(i<nxdin)
			w[0][i] 	= 0.0;
		else if(i>=nxdin&&i<(nx-nxdout))
			w[0][i] 	= 2.0*(psy0[0][i] 		- psy0[1][i])/(dx*dx); // top
		else
			w[0][i] 	= w[1][i];
		
		w[ny-1][i] 	= 2.0*(psy0[ny-1][i] 	- psy0[ny-2][i])/(dx*dx);	// bottom
	}
} // end of whichModel = 3.



if(whichModel==4){
// model 4 b.c.
	for(j=0;j<ny;j++){
		w[j][0] 		= 2.0*(psy0[j][0] 		- psy0[j][1])/(dx*dx); // left.
		w[j][nx-1] 	= w[j][nx-2]; // right.
	}
	for(i=0;i<nx;i++){
		if(i<nxdin)
			w[0][i] 	= 0.0;
		else 
			w[0][i] 	= 2.0*(psy0[0][i] 		- psy0[1][i])/(dx*dx); // top
		
		if(i<(nx-nxdin))
		w[ny-1][i] 	= 2.0*(psy0[ny-1][i] 	- psy0[ny-2][i])/(dx*dx);	// bottom
		else
		w[ny-1][i] = w[ny-2][i];
	}
} // end of whichModel = 4.




} // end of vorticityTransport.


// streamfunction.
double streamFunction(double **psy, double **psy0, double **w, int whichModel){

int i, j, it, tempd, tempb;
double residual;

for(it=0;it<1000;it++){
residual = 0.0;
for(i=1;i<nx-1;i++)
for(j=1;j<ny-1;j++){
psy[j][i] = betasor*0.25*(psy[j+1][i]+psy[j-1][i]+psy[j][i+1]+psy[j][i-1]+dx*dx*w[j][i])+(1.0 - betasor)*psy[j][i];
residual = residual + fabs(psy[j][i] - psy0[j][i]);
} // end of i,j loops

// blood
if(whichModel==-1){
	for(j=0;j<ny;j++){
		psy[j][nx-1] = ubr*Ly*((double)j/(double)(ny-1)); // right
		psy[j][0] 	= psy[j][1]; // left.
	}
	for(i=1;i<nx-1;i++){
		  psy[0][i] 	= 0.0; // ubr*Ly;
		  psy[ny-1][i] 	= ubr*Ly;
	}
} // end of whichModel = -1, blood.

if(whichModel==1){
	for(j=0;j<ny;j++){
	psy[j][nx-1] = psy[j][nx-2]; // right
	psy[j][0] 	= udl*Ly*((double)j/(double)(ny-1)); // left.
	}
	for(i=1;i<nx-1;i++){
		  psy[0][i] 	= 0.0;	  
		  psy[ny-1][i] 	= udl*Ly;
	}
} // end of whichModel = 1.

if(whichModel==2){
	for(j=0;j<ny;j++){
	psy[j][nx-1] 	= psy[j][nx-2]; // right
	psy[j][0] 		= udl*dx*(double)(nxdin-1); // left.
	}
	for(i=0;i<nx;i++){
	if(i<nxdin)
		psy[0][i] 		= udl*dx*(double)(nxdin-1) - udl*dx*(double)(i);
	else
		  psy[0][i] 	= 0.0;

		  psy[ny-1][i] 	= udl*dx*(double)(nxdin-1);
	}
} // end of whichModel = 2.


if(whichModel==3){
	for(j=0;j<ny;j++){
	psy[j][nx-1] 	= udl*dx*(double)(nxdin-1); // right
	psy[j][0] 		= udl*dx*(double)(nxdin-1); // left.
	}
	for(i=0;i<nx;i++){
	if(i<nxdin)
		psy[0][i] 		= udl*dx*(double)(nxdin-1) - udl*dx*(double)(i);
	else if(i>=nxdin&&i<(nx-nxdout))
		  psy[0][i] 	= 0.0;
	else
		  psy[0][i] = psy[1][i];

		  psy[ny-1][i] 	= udl*dx*(double)(nxdin-1);
	}
} // end of whichModel = 3.


if(whichModel==4){
	for(j=0;j<ny;j++){
	psy[j][nx-1] 	= 0.0; // right
	psy[j][0] 		= udl*dx*(double)(nxdin-1); // left.
	}
	for(i=0;i<nx;i++){
	if(i<nxdin)
		psy[0][i] 		= udl*dx*(double)(nxdin-1) - udl*dx*(double)(i);
	else if(i>=nxdin&&i<(nx-nxdout))
		  psy[0][i] 	= 0.0;

	if(i<(nx-nxdin))
		  psy[ny-1][i] 	= udl*dx*(double)(nxdin-1); // top.
	else
		psy[ny-1][i] = psy[ny-2][i];
	}
} // end of whichModel = 4.



// remplacer chaque fois le psyd0.
for(i=0;i<nx;i++) for(j=0;j<ny;j++) psy0[j][i] = psy[j][i];

if(residual<restol) break;
} // end of it loop.


return residual;
} // end of streamFunction

void velocities(double **psy, double **u, double **v, int whichModel){

int i, j;

for(j=1;j<ny-1;j++)
for(i=1;i<nx-1;i++){
	u[j][i] =  (psy[j+1][i] - psy[j-1][i])/(2.0*dy);
	v[j][i] = -(psy[j][i+1] - psy[j][i-1])/(2.0*dx);
}

if(whichModel==-1){
	for(j=0;j<ny;j++){
		u[j][0] 		= u[j][1];
		u[j][nx-1] 		= ubr; 
		v[j][0] 		= 0.0;
		v[j][nx-1] 		= v[j][nx-2];
	}
	for(i=1;i<nx-1;i++){
		u[0][i] 		= 0.0;
		u[ny-1][i] 		= 0.0;
		v[0][i] 		= 0.0;
		v[ny-1][i] 		= 0.0;
	} 
} // end of whichModel = 1.

if(whichModel==1){
	for(j=0;j<ny;j++){
		u[j][0] 		= udl;
		u[j][nx-1] 		= u[j][nx-2];
		v[j][0] 		= 0.0;
		v[j][nx-1] 		= v[j][nx-2];
	}
	for(i=1;i<nx-1;i++){
		u[0][i] 		= 0.0;
		u[ny-1][i] 		= 0.0;
		v[0][i] 		= 0.0;
		v[ny-1][i] 		= 0.0;
	} 
} // end of whichModel = 1.


if(whichModel==2){
	for(j=0;j<ny;j++){
		u[j][0] 		= 0.0;
		u[j][nx-1] 		= u[j][nx-2];
		v[j][0] 		= 0.0;
		v[j][nx-1] 		= v[j][nx-2];
	}
	for(i=0;i<nx;i++){
		u[0][i] 		= 0.0;
		u[ny-1][i] 		= 0.0;
		
		if(i<nxdin)
			v[0][i] 	= udl;
		else
			v[0][i] 	= 0.0;
			
		v[ny-1][i] 		= 0.0;
	} 
} // end of whichModel = 2.


if(whichModel==3){
	for(j=0;j<ny;j++){
		u[j][0] 		= 0.0;
		u[j][nx-1] 		= 0.0;
		v[j][0] 		= 0.0;
		v[j][nx-1] 		= 0.0;
	}
	for(i=0;i<nx;i++){
		u[0][i] 		= 0.0;
		u[ny-1][i] 		= 0.0;
		
		if(i<nxdin)
			v[0][i] 	= udl;
		else if(i>=nxdin&&i<(nx-nxdout))
			v[0][i] 	= 0.0;
		else
			v[0][i] = v[1][i];
			
		v[ny-1][i] 		= 0.0;
	} 
} // end of whichModel = 3.



if(whichModel==4){
	for(j=0;j<ny;j++){
		u[j][0] 		= 0.0;
		u[j][nx-1] 		= 0.0;
		v[j][0] 		= 0.0;
		v[j][nx-1] 		= 0.0;
	}
	for(i=0;i<nx;i++){
		u[0][i] 		= 0.0;
		u[ny-1][i] 		= 0.0;
		
		if(i<nxdin)
			v[0][i] 	= udl;
		else
			v[0][i] 	= 0.0;

		if(i<(nx-nxdout))		
		v[ny-1][i] 		= 0.0;
		else
		v[ny-1][i] = v[ny-2][i];
	} 
} // end of whichModel = 4.




return;
}


double pressure(double **p, double **psy, double **u, double **v, int whichModel){

int i, j, itt;
double residual, diff, d2psydx2, d2psydy2, prod, prod2, d2psydxdy, tot, sor_adjustment;
double rho, mu;
if(whichModel==-1){
rho = rhob; mu = mub;
}
else{
rho = rhod; mu = mud;
}


for(itt=0;itt<1000000000;itt++){

residual = 0.0;
for(j=1;j<ny-1;j++)
for(i=1;i<nx-1;i++){

diff 			= 0.25*(p[j][i+1] + p[j][i-1] + p[j+1][i] + p[j-1][i]);
d2psydx2 	=          (psy[j][i+1] + psy[j][i-1] - 2.0*psy[j][i]);
d2psydy2 	=         (psy[j+1][i] + psy[j-1][i] - 2.0*psy[j][i]);
prod 		= 0.25*d2psydx2*d2psydy2;
d2psydxdy 	=          (psy[j][i+1] + psy[j][i-1] + psy[j+1][i] + psy[j-1][i] - 4.0*psy[j][i]);
prod2 		= 0.25*d2psydxdy*d2psydxdy;

tot 			= diff - 2.0*rho*(prod - prod2)/(dx*dx);
sor_adjustment = alpha * (tot - p[j][i]);

p[j][i] 		= sor_adjustment  + tot; // because dx = dy.

residual 		= residual + fabs(p[j][i] - tot);

} // end of i, j.


if(whichModel==-1){
	// blood b.c.
	for(j=0;j<ny;j++)		p[j][0] = pbL; // left.Ding, porosity paper table 1 rectified.
	for(i=0;i<nx;i++){
		p[0][i] 			= p[1][i]; // top, solid boundary.
		p[ny-1][i] 			= p[ny-2][i]; // bottom, solid boundary.
	}
	for(j=0;j<ny;j++)		p[j][nx-1] = p[j][nx-2] - dx*mub*u[j][nx-1]/kbx;
} // end of whichModel.



if(whichModel==1){
	// dialysate b.c.
	for(j=0;j<ny;j++)		p[j][0] = p[j][1] + dx*mud*u[j][0]/kdx; // left.Ding, porosity paper table 1 rectified.
	for(i=1;i<nx-1;i++){
		p[0][i] 			= p[1][i]; 		// top, solid boundary.
		p[ny-1][i] 			= p[ny-2][i]; 	// bottom, solid boundary.
	}
	for(j=0;j<ny;j++)		p[j][nx-1] = pdR; 	// right, Dirichlet.
} // end of whichModel 1.

if(whichModel==2){
	// dialysate b.c.
	for(j=0;j<ny;j++)		p[j][0] = p[j][1]; // left.
	for(i=0;i<nx;i++){
	if(i<nxdin)
		p[0][i] 			= p[1][i] + dx*mud*u[0][j]/kdy; 
	else
		p[0][i] 			= p[1][i]; 			// top, solid boundary.
		
		p[ny-1][i] 			= p[ny-2][i]; 		// bottom, solid boundary.
	}
	for(j=0;j<ny;j++)		p[j][nx-1] = pdR; 	// right, Dirichlet.
} // end of whichModel 2.

if(whichModel==3){
	// dialysate b.c.
	for(j=0;j<ny;j++)		p[j][0] = p[j][1]; // left.
	for(i=0;i<nx;i++){
	if(i<nxdin)
		p[0][i] 			= p[1][i] + dx*mud*v[0][j]/kdy; 
	else if(i>nxdin&&i<(nx-nxdout))
		p[0][i] 			= p[1][i]; 			// top, solid boundary.
	else
		p[0][i] = pdR;
		
		p[ny-1][i] 			= p[ny-2][i]; 		// bottom, solid boundary.
	}
	for(j=0;j<ny;j++)		p[j][nx-1] = p[j][nx-2]; 	// right.
} // end of whichModel 3.



if(whichModel==4){
	// dialysate b.c.
	for(j=0;j<ny;j++)		p[j][0] = p[j][1]; // left.
	for(i=0;i<nx;i++){
	if(i<nxdin)
		p[0][i] 			= p[1][i] + dx*mud*v[0][j]/kdy; 
	else 
		p[0][i] 			= p[1][i]; 			// top, solid boundary.

	if(i<(nx-nxdout))	
		p[ny-1][i] 			= p[ny-2][i]; 		// bottom, solid boundary.
	else
		p[ny-1][i] = pdR;
	}
	for(j=0;j<ny;j++)		p[j][nx-1] = p[j][nx-2]; 	// right.
} // end of whichModel 4.




if(itt%1000==0) printf("pressure residual %f for model %d\n", residual, whichModel);
// for model 1, restol/50 was enough.
if(residual<restol/50.0&&itt>100) break;

} // end of itt.

return residual;
} // end of pressure.


void ureabcs(double **c, int whichModel){

int i, j;

if(whichModel==-1){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= cbru;
		c[j][0] 	= c[j][1]; // this must be neumann when you start the dialysate-blood transport.
	}
	for(i=0;i<nx;i++){
		c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = -1.

if(whichModel==1){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= c[j][nx-2];
		c[j][0] 	= cdlu;
	}
	for(i=0;i<nx;i++){
		c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = 1.


if(whichModel==2||whichModel==3||whichModel==4){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= c[j][nx-2];
		c[j][0] 	= c[j][1];
	}
	for(i=0;i<nx;i++){
		if(i<nxdin)
			c[0][i] = cdlu;
		else
			c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = 2,3,4.




} // end of urea b.c.s


/*
The Lax-Wendroff on its own will not work.
First, split the advection and diffusion. Do advection first.
In the advection, do the x direction once, then do the y direction. see:
https://scicomp.stackexchange.com/questions/28645/finite-volume-piecewise-linear-2d-advection-develops-instability
*/
double urea(double **c, double **c0, double **u, double **v, double diffx, double diffy, int whichModel){

int i, j, iter;
double **cp1, **cp2;
double residual, conv1, conv2, diffterm;

cp1 = doubleAllocate( ny, nx );
cp2 = doubleAllocate( ny, nx );

for(i=0;i<nx;i++) for(j=0;j<ny;j++){ cp1[j][i] = cp2[j][i] = 0.0; }

for(iter=0;iter<500000;iter++){
residual = 0.0;

// x sweep.
for(j=0;j<ny;j++)
for(i=1;i<nx-1;i++){
	conv1 	= DELTAT*			( -u[j][i]*(c0[j][i+1] - c0[j][i-1])/(2.0*dx) );
	conv2 	= DELTAT*DELTAT*	(u[j][i]*u[j][i]*(c0[j][i+1]+c0[j][i-1]-2.0*c0[j][i])/(2.0*dx*dx));
	diffterm 	= DELTAT*(diffx*(c0[j][i+1]+c0[j][i-1]-2.0*c0[j][i])/(dx*dx)); 
	cp1[j][i]	= c0[j][i] + conv1; // + conv2; // + diffterm*betasor + (1.0 - betasor)*c0[j][i];
} // end of x sweep.
ureabcs(cp1, whichModel);

// y sweep.
for(i=0;i<nx;i++)
for(j=1;j<ny-1;j++){
	conv1 	= DELTAT*			( -v[j][i]*(cp1[j+1][i] - cp1[j-1][i])/(2.0*dx) );
	conv2 	= DELTAT*DELTAT*	(v[j][i]*v[j][i]*(cp1[j+1][i]+cp1[j-1][i]-2.0*cp1[j][i])/(2.0*dx*dx));
	diffterm 	= DELTAT*(diffy*(cp1[j+1][i]+cp1[j-1][i]-2.0*cp1[j][i])/(dx*dx));
	cp2[j][i]	= cp1[j][i] + conv1; // + conv2; // + diffterm*betasor + (1.0 - betasor)*cp1[j][i];
} // end of y sweep.
ureabcs(cp2, whichModel);

// do diffusion with SOR.
for(i=1;i<nx-1;i++)
for(j=1;j<ny-1;j++){
	diffterm 	= DELTAT*(diffx*(cp2[j][i+1]+cp2[j][i-1]-2.0*cp2[j][i])/(dx*dx)+diffy*(cp2[j+1][i]+cp2[j-1][i]-2.0*cp2[j][i])/(dx*dx));
	c[j][i] 	= cp2[j][i] + betasor*diffterm + (1.0 - betasor)*c[j][i];
//	c[j][i] 	= cp2[j][i] + diffterm;
} // end of convection loops.


ureabcs(c, whichModel);


for(j=0;j<ny;j++) for(i=0;i<nx;i++){
 residual = residual + fabs(c[j][i] - c0[j][i]);
 c0[j][i] = c[j][i];
}

printf("urea residual: %20.100f in model %d %d\n", residual, whichModel, iter);
 if(residual<1.0e-20) break;
} // end of iter.

deallocate(cp1, ny, nx);
deallocate(cp2, ny, nx);

return residual;
} // end of urea.


double ureaOLD(double **c, double **c0, double **u, double **v, double diffx, double diffy, int whichModel){
int i, j, iter;
double conv1, conv2, conv3, conv4, conv5, conv6, conv7, diffterm;
double cp[ny][nx];
double residual;

for(iter=0;iter<50000;iter++){
residual = 0.0;
for(i=1;i<nx-1;i++)
for(j=1;j<ny-1;j++){
	conv1 	= DELTAT*			( -u[j][i]*(c0[j][i+1] - c0[j][i-1])/(2.0*dx) - v[j][i]*(c0[j+1][i] - c0[j-1][i])/(2.0*dx) );
	conv2 	= DELTAT*DELTAT*	(u[j][i]*u[j][i]*(c0[j][i+1]+c0[j][i-1]-2.0*c0[j][i])/(dx*dx));
	conv3 	= DELTAT*DELTAT*	(v[j][i]*v[j][i]*(c0[j+1][i]+c0[j-1][i]-2.0*c0[j][i])/(dx*dx));
	conv4 	= DELTAT*DELTAT*	2.0*u[j][i]*v[j][i]*(c0[j+1][i+1]+c0[j-1][i-1]-c0[j-1][i+1]-c0[j+1][i-1])/(4.0*dx*dx);
	diffterm 	= DELTAT*			(diffx*(c0[j][i+1]+c0[j][i-1]-2.0*c0[j][i])/(dx*dx)+diffy*(c0[j+1][i]+c0[j-1][i]-2.0*c0[j][i])/(dx*dx));
	c[j][i] 	= c0[j][i] + conv1 + conv2 + conv3 + conv4 + betasor*diffterm + (1.0 - betasor)*c0[j][i];
//	c[j][i] 	= c0[j][i] + conv1 + diffterm;
} // end of convection loops.

if(whichModel==-1){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= cbru;
		c[j][0] 	= c[j][1]; // this must be neumann when you start the dialysate-blood transport.
	}
	for(i=0;i<nx;i++){
		c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = -1.

if(whichModel==1){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= c[j][nx-2];
		c[j][0] 	= cdlu;
	}
	for(i=0;i<nx;i++){
		c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = 1.


if(whichModel==2||whichModel==3||whichModel==4){
	for(j=0;j<ny;j++){
		c[j][nx-1] 	= c[j][nx-2];
		c[j][0] 	= c[j][1];
	}
	for(i=0;i<nx;i++){
		if(i<nxdin)
			c[0][i] = cdlu;
		else
			c[0][i] 	= c[1][i];
		c[ny-1][i] 	= c[ny-2][i];
	}
} // end of whichModel = 2,3,4.




for(j=0;j<ny;j++) for(i=0;i<nx;i++){
 residual = residual + fabs(c[j][i] - c0[j][i]);
 c0[j][i] = c[j][i];
}

printf("urea residual: %20.100f in model %d %d\n", residual, whichModel, iter);
// sleep(1);
 if(residual<1.0e-15) break;
} // end of iter loops.

return residual;
} // end of ureaOLD.

void JvS1(double **Jv, double **S1, double **pb, double **pd){
int i, j;
for(j=0;j<ny;j++)
for(i=0;i<nx;i++){
	Jv[j][i] 				= Lp*(pb[j][i] - pd[j][i]);
	if(Jv[j][i]>0.0)	S1[j][i] 	= betaa*Jv[j][i]*rhob;
	else			S1[j][i] 	=  betaa*Jv[j][i]*rhod;
} // end of i, j loops.	
} // end of JvS1 function.


