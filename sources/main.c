/*
Sanjay R Kharche.
28 April 2022. I tried working things on MATLAB.
The matlab tests are taking more than 4 hours. It is now worth doing the C programming for a few days.
Serial program.

see example in matlab:
https://sehenaoga.github.io/projects/NS_vorticity.html
*/

#include "srkDialyser.h"

int main(int argc, char *argv[]){

FILE *fpsyd, *fwd, *fud, *fvd, *fpd;
FILE *fpsyb, *fwb, *fub, *fvb, *fpb;
FILE *fcdu, *fcbu;
FILE *fJv, *fS1;

double **psyd, **psyd0, **wd, **wd0, **ud, **vd, **pd;
double **psyb, **psyb0, **wb, **wb0, **ub, **vb, **pb;
double **cdu, **cdu0, **cbu, **cbu0;
double **Jv, **S1;

int i, j, t, t_vpsi, whichModel;
double reswd, respsyd, respd;
double reswb, respsyb, respb;
double rescbu;

/* whichModel 
-1: blood.
0: no such case.
1: straight pipe left to right, dialysate.
2: i = 0 to nxdin inlet, outlet at right.
3:
4:
5:
*/
whichModel = 1; 

psyd 	= doubleAllocate( ny, nx );
psyd0 	= doubleAllocate( ny, nx );
wd 		= doubleAllocate( ny, nx );
wd0 	= doubleAllocate( ny, nx );
ud 		= doubleAllocate( ny, nx );
vd	 	= doubleAllocate( ny, nx );
pd	 	= doubleAllocate( ny, nx );

psyb 	= doubleAllocate( ny, nx );
psyb0 	= doubleAllocate( ny, nx );
wb 		= doubleAllocate( ny, nx );
wb0 	= doubleAllocate( ny, nx );
ub 		= doubleAllocate( ny, nx );
vb	 	= doubleAllocate( ny, nx );
pb	 	= doubleAllocate( ny, nx );

cbu	 	= doubleAllocate( ny, nx );
cbu0 	= doubleAllocate( ny, nx );
cdu 		= doubleAllocate( ny, nx );
cdu0 	= doubleAllocate( ny, nx );

Jv	 	= doubleAllocate( ny, nx );
S1	 	= doubleAllocate( ny, nx );


// initialize(psyd, psyd0); // may not be needed, reduces editing by not calling this function.

for(t_vpsi = 0; t_vpsi<100000; t_vpsi++){
	vorticityTransport(wd, wd0, psyd0, S1, whichModel, 0.0);
	vorticityTransport(wb, wb0, psyb0, S1, -1, 0.0);
	reswd 	= residuals(wd, wd0, 1);
	reswb 	= residuals(wb, wb0, 1);
	respsyd 	= streamFunction(psyd, psyd0, wd, whichModel);
	respsyb 	= streamFunction(psyb, psyb0, wb, -1);
	printf("%f %f %f %f\n", reswd, respsyd, reswb, respsyb);
	if(reswd<restol&&respsyd<restol&&t_vpsi>1000) break;
} // end of t loop.
//*************************************************************************
//*************************************************************************

// velocity (one iterate) after psi and w have converged.
velocities(psyd, ud, vd, whichModel);
velocities(psyb, ub, vb, -1);
// urea (many iterations). 2 May 2022. Solving the advection diffusion is difficult. I will do other things first.
// rescbu = urea(cbu, cbu0, ub, vb, cbdiffux, cbdiffuy, -1);
// rescbu = urea(cdu, cdu0, ud, vd, cddiffux, cddiffuy, whichModel);
// pressure (many iterations) after psi and w have converged.
respb = pressure(pb, psyb, ub, vb, -1);
respd = pressure(pd, psyd, ud, vd, whichModel);

// get the Jv and S1.
JvS1(Jv, S1, pb, pd);


// the psy and w may be settled now. Do the darcy.
for(t_vpsi = 0; t_vpsi<10; t_vpsi++){
	printf("Now you can do inhomogeneous terms, which also makes the pdes unstable. after conditioning.\n");
}
	
// outputs.
fpsyd 	= fopen("psyd.dat","w");
fwd 		= fopen("wd.dat","w");
fud 		= fopen("ud.dat","w");
fpd 		= fopen("pd.dat","w");
// fcdu 	= fopen("cdu.dat","w");

fpsyb 	= fopen("psyb.dat","w");
fwb 		= fopen("wb.dat","w");
fub 		= fopen("ub.dat","w");
fpb 		= fopen("pb.dat","w");
// fcbu 	= fopen("cbu.dat","w");

fJv = fopen("Jv.dat","w");
fS1 = fopen("S1.dat","w");

for(j=0;j<ny;j++)
for(i=0;i<nx;i++){
	fprintf(fpsyd, 	"%d %d %20.50f\n", 			i, j, psyd[j][i])	;
	fprintf(fwd, 	"%d %d %20.50f\n", 			i, j, wd[j][i])	;
	fprintf(fud, 	"%d %d %20.50f %20.50f\n", 	i, j, ud[j][i], vd[j][i] )	;	
	fprintf(fpd, 	"%d %d %20.50f\n", 			i, j, pd[j][i]/1333.3)	;
//	fprintf(fcdu, 	"%d %d %20.50f\n", 			i, j, cdu[j][i])	;	
	
	fprintf(fpsyb, 	"%d %d %20.50f\n", 			i, j, psyb[j][i])	;
	fprintf(fwb, 	"%d %d %20.50f\n", 			i, j, wb[j][i])	;
	fprintf(fub, 	"%d %d %20.50f %20.50f\n", 	i, j, ub[j][i], vb[j][i] )	;	
	fprintf(fpb, 	"%d %d %20.50f\n", 			i, j, pb[j][i]/1333.3)	;
//	fprintf(fcbu, 	"%d %d %20.50f\n", 			i, j, cbu[j][i])	;

fprintf(fJv, "%d %d %20.50f\n", 			i, j, Jv[j][i])	;
fprintf(fS1, "%d %d %20.50f\n", 			i, j, S1[j][i])	;

}

fclose(fpsyd);
fclose(fwd);
fclose(fud);
fclose(fpd);
// fclose(fcdu);

fclose(fpsyb);
fclose(fwb);
fclose(fub);
fclose(fpb);
// fclose(fcbu);

fclose(fJv);
fclose(fS1);


return 0;
} // end of main.
