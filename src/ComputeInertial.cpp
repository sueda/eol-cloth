#include <cmath>
void ComputeInertial(
   const double *xa, // [input 3x1] World position of vertex A
   const double *xb, // [input 3x1] World position of vertex B
   const double *xc, // [input 3x1] World position of vertex C
   const double *Xa, // [input 2x1] Material position of vertex A
   const double *Xb, // [input 2x1] Material position of vertex B
   const double *Xc, // [input 2x1] Material position of vertex C
   const double *g,  // [input 3x1] 3D gravity vector
   double rho,       // [input 1x1] Density (mass per area)
   double *W,        // [output 1x1] Gravitational potential energy
   double *f,        // [output 9x1] Gravity force vector
   double *M)        // [output 9x9] Inertia matrix
{
double xax = xa[0];
double xay = xa[1];
double xaz = xa[2];
double xbx = xb[0];
double xby = xb[1];
double xbz = xb[2];
double xcx = xc[0];
double xcy = xc[1];
double xcz = xc[2];
double Xax = Xa[0];
double Xay = Xa[1];
double Xbx = Xb[0];
double Xby = Xb[1];
double Xcx = Xc[0];
double Xcy = Xc[1];
double gx = g[0];
double gy = g[1];
double gz = g[2];
double t8 = rho * (Xax * Xby - Xax * Xcy - Xbx * Xay + Xcx * Xay + Xbx * Xcy - Xcx * Xby);
double W00 = -t8 * (gx * (xbx - xcx) / 6 + gy * (xby - xcy) / 6 + gz * (xbz - xcz) / 6 + gx * (xax - xcx) / 6 + gy * (xay - xcy) / 6 + gz * (xaz - xcz) / 6 + gx * xcx / 2 + gy * xcy / 2 + gz * xcz / 2);
double f01 = t8 * gx / 6;
double f02 = t8 * gy / 6;
double f03 = t8 * gz / 6;
double f04 = f01;
double f05 = f02;
double f06 = f03;
double f07 = f04;
double f08 = f05;
double f09 = f06;
double M0101 = t8 / 12;
double M0102 = 0;
double M0103 = 0;
double M0104 = t8 / 24;
double M0105 = 0;
double M0106 = 0;
double M0107 = M0104;
double M0108 = 0;
double M0109 = 0;
double M0201 = 0;
double M0202 = M0101;
double M0203 = 0;
double M0204 = 0;
double M0205 = M0107;
double M0206 = 0;
double M0207 = 0;
double M0208 = M0205;
double M0209 = 0;
double M0301 = 0;
double M0302 = 0;
double M0303 = M0202;
double M0304 = 0;
double M0305 = 0;
double M0306 = M0208;
double M0307 = 0;
double M0308 = 0;
double M0309 = M0306;
double M0401 = M0309;
double M0402 = 0;
double M0403 = 0;
double M0404 = M0303;
double M0405 = 0;
double M0406 = 0;
double M0407 = M0401;
double M0408 = 0;
double M0409 = 0;
double M0501 = 0;
double M0502 = M0407;
double M0503 = 0;
double M0504 = 0;
double M0505 = M0404;
double M0506 = 0;
double M0507 = 0;
double M0508 = M0502;
double M0509 = 0;
double M0601 = 0;
double M0602 = 0;
double M0603 = M0508;
double M0604 = 0;
double M0605 = 0;
double M0606 = M0505;
double M0607 = 0;
double M0608 = 0;
double M0609 = M0603;
double M0701 = M0609;
double M0702 = 0;
double M0703 = 0;
double M0704 = M0701;
double M0705 = 0;
double M0706 = 0;
double M0707 = M0606;
double M0708 = 0;
double M0709 = 0;
double M0801 = 0;
double M0802 = M0704;
double M0803 = 0;
double M0804 = 0;
double M0805 = M0802;
double M0806 = 0;
double M0807 = 0;
double M0808 = M0707;
double M0809 = 0;
double M0901 = 0;
double M0902 = 0;
double M0903 = M0805;
double M0904 = 0;
double M0905 = 0;
double M0906 = M0903;
double M0907 = 0;
double M0908 = 0;
double M0909 = M0808;
W[0]=W00;
f[0]=f01; f[1]=f02; f[2]=f03; f[3]=f04; f[4]=f05; f[5]=f06; f[6]=f07; f[7]=f08; f[8]=f09; 
M[ 0*9+ 0]=M0101; M[ 0*9+ 1]=M0102; M[ 0*9+ 2]=M0103; M[ 0*9+ 3]=M0104; M[ 0*9+ 4]=M0105; M[ 0*9+ 5]=M0106; M[ 0*9+ 6]=M0107; M[ 0*9+ 7]=M0108; M[ 0*9+ 8]=M0109; 
M[ 1*9+ 0]=M0201; M[ 1*9+ 1]=M0202; M[ 1*9+ 2]=M0203; M[ 1*9+ 3]=M0204; M[ 1*9+ 4]=M0205; M[ 1*9+ 5]=M0206; M[ 1*9+ 6]=M0207; M[ 1*9+ 7]=M0208; M[ 1*9+ 8]=M0209; 
M[ 2*9+ 0]=M0301; M[ 2*9+ 1]=M0302; M[ 2*9+ 2]=M0303; M[ 2*9+ 3]=M0304; M[ 2*9+ 4]=M0305; M[ 2*9+ 5]=M0306; M[ 2*9+ 6]=M0307; M[ 2*9+ 7]=M0308; M[ 2*9+ 8]=M0309; 
M[ 3*9+ 0]=M0401; M[ 3*9+ 1]=M0402; M[ 3*9+ 2]=M0403; M[ 3*9+ 3]=M0404; M[ 3*9+ 4]=M0405; M[ 3*9+ 5]=M0406; M[ 3*9+ 6]=M0407; M[ 3*9+ 7]=M0408; M[ 3*9+ 8]=M0409; 
M[ 4*9+ 0]=M0501; M[ 4*9+ 1]=M0502; M[ 4*9+ 2]=M0503; M[ 4*9+ 3]=M0504; M[ 4*9+ 4]=M0505; M[ 4*9+ 5]=M0506; M[ 4*9+ 6]=M0507; M[ 4*9+ 7]=M0508; M[ 4*9+ 8]=M0509; 
M[ 5*9+ 0]=M0601; M[ 5*9+ 1]=M0602; M[ 5*9+ 2]=M0603; M[ 5*9+ 3]=M0604; M[ 5*9+ 4]=M0605; M[ 5*9+ 5]=M0606; M[ 5*9+ 6]=M0607; M[ 5*9+ 7]=M0608; M[ 5*9+ 8]=M0609; 
M[ 6*9+ 0]=M0701; M[ 6*9+ 1]=M0702; M[ 6*9+ 2]=M0703; M[ 6*9+ 3]=M0704; M[ 6*9+ 4]=M0705; M[ 6*9+ 5]=M0706; M[ 6*9+ 6]=M0707; M[ 6*9+ 7]=M0708; M[ 6*9+ 8]=M0709; 
M[ 7*9+ 0]=M0801; M[ 7*9+ 1]=M0802; M[ 7*9+ 2]=M0803; M[ 7*9+ 3]=M0804; M[ 7*9+ 4]=M0805; M[ 7*9+ 5]=M0806; M[ 7*9+ 6]=M0807; M[ 7*9+ 7]=M0808; M[ 7*9+ 8]=M0809; 
M[ 8*9+ 0]=M0901; M[ 8*9+ 1]=M0902; M[ 8*9+ 2]=M0903; M[ 8*9+ 3]=M0904; M[ 8*9+ 4]=M0905; M[ 8*9+ 5]=M0906; M[ 8*9+ 6]=M0907; M[ 8*9+ 7]=M0908; M[ 8*9+ 8]=M0909; 
}
