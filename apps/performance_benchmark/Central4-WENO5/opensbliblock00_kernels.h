#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
 void opensbliblock00Kernel037(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
ACC<double> &rhou2_B0, ACC<double> &x0_B0, ACC<double> &x1_B0, ACC<double> &x2_B0, const int *idx)
{
   double p = 0.0;
   double r = 0.0;
   double u0 = 0.0;
   double u1 = 0.0;
   double u2 = 0.0;
   double x0 = 0.0;
   double x1 = 0.0;
   double x2 = 0.0;
   x0 = Delta0block0*idx[0];

   x1 = Delta1block0*idx[1];

   x2 = Delta2block0*idx[2];

   u0 = cos(x1)*cos(x2)*sin(x0);

   u1 = -cos(x0)*cos(x2)*sin(x1);

   u2 = 0.0;

   p = (2.0 + cos(2.0*x2))*(0.0625*cos(2.0*x0) + 0.0625*cos(2.0*x1)) + 1.0/((Minf*Minf)*gama);

   r = gama*p*(Minf*Minf);

   rho_B0(0,0,0) = r;

   rhou0_B0(0,0,0) = r*u0;

   rhou1_B0(0,0,0) = r*u1;

   rhou2_B0(0,0,0) = r*u2;

   rhoE_B0(0,0,0) = p/(-1 + gama) + 0.5*r*((u0*u0) + (u1*u1) + (u2*u2));

   x0_B0(0,0,0) = x0;

   x1_B0(0,0,0) = x1;

   x2_B0(0,0,0) = x2;

}

void opensbliblock00Kernel005(const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, ACC<double> &u0_B0)
{
   u0_B0(0,0,0) = rhou0_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel007(const ACC<double> &rho_B0, const ACC<double> &rhou1_B0, ACC<double> &u1_B0)
{
   u1_B0(0,0,0) = rhou1_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel009(const ACC<double> &rho_B0, const ACC<double> &rhou2_B0, ACC<double> &u2_B0)
{
   u2_B0(0,0,0) = rhou2_B0(0,0,0)/rho_B0(0,0,0);

}

 void opensbliblock00Kernel019(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &p_B0)
{
    p_B0(0,0,0) = (-1 + gama)*(-(1.0/2.0)*(u0_B0(0,0,0)*u0_B0(0,0,0))*rho_B0(0,0,0) -
      (1.0/2.0)*(u1_B0(0,0,0)*u1_B0(0,0,0))*rho_B0(0,0,0) - (1.0/2.0)*(u2_B0(0,0,0)*u2_B0(0,0,0))*rho_B0(0,0,0) +
      rhoE_B0(0,0,0));

}

 void opensbliblock00Kernel016(const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const ACC<double> &rho_B0,
ACC<double> &H_B0)
{
   H_B0(0,0,0) = (p_B0(0,0,0) + rhoE_B0(0,0,0))/rho_B0(0,0,0);

}

void opensbliblock00Kernel022(const ACC<double> &p_B0, const ACC<double> &rho_B0, ACC<double> &T_B0)
{
   T_B0(0,0,0) = (Minf*Minf)*gama*p_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel023(const ACC<double> &T_B0, ACC<double> &mu_B0)
{
   mu_B0(0,0,0) = 1.4042*T_B0(0,0,0)*sqrt(T_B0(0,0,0))/(0.40417 + T_B0(0,0,0));

}

void opensbliblock00Kernel004(const ACC<double> &u0_B0, ACC<double> &wk0_B0)
{
    wk0_B0(0,0,0) = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel006(const ACC<double> &u1_B0, ACC<double> &wk1_B0)
{
    wk1_B0(0,0,0) = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel008(const ACC<double> &u2_B0, ACC<double> &wk2_B0)
{
    wk2_B0(0,0,0) = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel010(const ACC<double> &u0_B0, ACC<double> &wk3_B0)
{
    wk3_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
      ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel011(const ACC<double> &u1_B0, ACC<double> &wk4_B0)
{
    wk4_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
      ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel012(const ACC<double> &u2_B0, ACC<double> &wk5_B0)
{
    wk5_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
      ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel013(const ACC<double> &u0_B0, ACC<double> &wk6_B0)
{
    wk6_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel014(const ACC<double> &u1_B0, ACC<double> &wk7_B0)
{
    wk7_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel015(const ACC<double> &u2_B0, ACC<double> &wk8_B0)
{
    wk8_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

}

 void opensbliblock00Kernel029(const ACC<double> &H_B0, const ACC<double> &p_B0, const ACC<double> &rho_B0, const
ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, const ACC<double> &wk0_B0, const ACC<double> &wk1_B0, const ACC<double>
&wk2_B0, const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, const ACC<double> &wk5_B0, const ACC<double> &wk6_B0,
const ACC<double> &wk7_B0, const ACC<double> &wk8_B0, ACC<double> &Residual0_B0, ACC<double> &Residual1_B0, ACC<double>
&Residual2_B0, ACC<double> &Residual3_B0, ACC<double> &Residual4_B0)
{
   double d1_H_dx = 0.0;
   double d1_H_dy = 0.0;
   double d1_H_dz = 0.0;
   double d1_Hrho_dx = 0.0;
   double d1_Hrho_dy = 0.0;
   double d1_Hrho_dz = 0.0;
   double d1_Hrhou0_dx = 0.0;
   double d1_Hrhou1_dy = 0.0;
   double d1_Hrhou2_dz = 0.0;
   double d1_Hu0_dx = 0.0;
   double d1_Hu1_dy = 0.0;
   double d1_Hu2_dz = 0.0;
   double d1_p_dx = 0.0;
   double d1_p_dy = 0.0;
   double d1_p_dz = 0.0;
   double d1_rho_dx = 0.0;
   double d1_rho_dy = 0.0;
   double d1_rho_dz = 0.0;
   double d1_rhou0_dx = 0.0;
   double d1_rhou0_dy = 0.0;
   double d1_rhou0_dz = 0.0;
   double d1_rhou0u0_dx = 0.0;
   double d1_rhou0u1_dx = 0.0;
   double d1_rhou0u2_dx = 0.0;
   double d1_rhou1_dx = 0.0;
   double d1_rhou1_dy = 0.0;
   double d1_rhou1_dz = 0.0;
   double d1_rhou1u0_dy = 0.0;
   double d1_rhou1u1_dy = 0.0;
   double d1_rhou1u2_dy = 0.0;
   double d1_rhou2_dx = 0.0;
   double d1_rhou2_dy = 0.0;
   double d1_rhou2_dz = 0.0;
   double d1_rhou2u0_dz = 0.0;
   double d1_rhou2u1_dz = 0.0;
   double d1_rhou2u2_dz = 0.0;
   double d1_u0u0_dx = 0.0;
   double d1_u0u1_dx = 0.0;
   double d1_u0u1_dy = 0.0;
   double d1_u0u2_dx = 0.0;
   double d1_u0u2_dz = 0.0;
   double d1_u1u1_dy = 0.0;
   double d1_u1u2_dy = 0.0;
   double d1_u1u2_dz = 0.0;
   double d1_u2u2_dz = 0.0;
    d1_H_dx = (-(2.0/3.0)*H_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0) + ((1.0/12.0))*H_B0(-2,0,0) +
      ((2.0/3.0))*H_B0(1,0,0))*invDelta0block0;

    d1_Hrho_dx = (-(2.0/3.0)*H_B0(-1,0,0)*rho_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*rho_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*rho_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*rho_B0(1,0,0))*invDelta0block0;

    d1_Hrhou0_dx = (-(2.0/3.0)*H_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_Hu0_dx = (-(2.0/3.0)*H_B0(-1,0,0)*u0_B0(-1,0,0) - (1.0/12.0)*H_B0(2,0,0)*u0_B0(2,0,0) +
      ((1.0/12.0))*H_B0(-2,0,0)*u0_B0(-2,0,0) + ((2.0/3.0))*H_B0(1,0,0)*u0_B0(1,0,0))*invDelta0block0;

    d1_p_dx = (-(2.0/3.0)*p_B0(-1,0,0) - (1.0/12.0)*p_B0(2,0,0) + ((1.0/12.0))*p_B0(-2,0,0) +
      ((2.0/3.0))*p_B0(1,0,0))*invDelta0block0;

    d1_rho_dx = (-(2.0/3.0)*rho_B0(-1,0,0) - (1.0/12.0)*rho_B0(2,0,0) + ((1.0/12.0))*rho_B0(-2,0,0) +
      ((2.0/3.0))*rho_B0(1,0,0))*invDelta0block0;

    d1_rhou0_dx = (-(2.0/3.0)*rhou0_B0(-1,0,0) - (1.0/12.0)*rhou0_B0(2,0,0) + ((1.0/12.0))*rhou0_B0(-2,0,0) +
      ((2.0/3.0))*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u1_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u1_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou0u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0)*rhou0_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0)*rhou0_B0(2,0,0) +
      ((1.0/12.0))*u2_B0(-2,0,0)*rhou0_B0(-2,0,0) + ((2.0/3.0))*u2_B0(1,0,0)*rhou0_B0(1,0,0))*invDelta0block0;

    d1_rhou1_dx = (-(2.0/3.0)*rhou1_B0(-1,0,0) - (1.0/12.0)*rhou1_B0(2,0,0) + ((1.0/12.0))*rhou1_B0(-2,0,0) +
      ((2.0/3.0))*rhou1_B0(1,0,0))*invDelta0block0;

    d1_rhou2_dx = (-(2.0/3.0)*rhou2_B0(-1,0,0) - (1.0/12.0)*rhou2_B0(2,0,0) + ((1.0/12.0))*rhou2_B0(-2,0,0) +
      ((2.0/3.0))*rhou2_B0(1,0,0))*invDelta0block0;

    d1_u0u0_dx = (-(2.0/3.0)*(u0_B0(-1,0,0)*u0_B0(-1,0,0)) - (1.0/12.0)*(u0_B0(2,0,0)*u0_B0(2,0,0)) +
      ((1.0/12.0))*(u0_B0(-2,0,0)*u0_B0(-2,0,0)) + ((2.0/3.0))*(u0_B0(1,0,0)*u0_B0(1,0,0)))*invDelta0block0;

    d1_u0u1_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*u1_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*u1_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*u1_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*u1_B0(1,0,0))*invDelta0block0;

    d1_u0u2_dx = (-(2.0/3.0)*u0_B0(-1,0,0)*u2_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0)*u2_B0(2,0,0) +
      ((1.0/12.0))*u0_B0(-2,0,0)*u2_B0(-2,0,0) + ((2.0/3.0))*u0_B0(1,0,0)*u2_B0(1,0,0))*invDelta0block0;

    d1_H_dy = (-(2.0/3.0)*H_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0) + ((1.0/12.0))*H_B0(0,-2,0) +
      ((2.0/3.0))*H_B0(0,1,0))*invDelta1block0;

    d1_Hrho_dy = (-(2.0/3.0)*H_B0(0,-1,0)*rho_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*rho_B0(0,2,0) +
      ((1.0/12.0))*H_B0(0,-2,0)*rho_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*rho_B0(0,1,0))*invDelta1block0;

    d1_Hrhou1_dy = (-(2.0/3.0)*H_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*rhou1_B0(0,2,0) +
      ((1.0/12.0))*H_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

    d1_Hu1_dy = (-(2.0/3.0)*H_B0(0,-1,0)*u1_B0(0,-1,0) - (1.0/12.0)*H_B0(0,2,0)*u1_B0(0,2,0) +
      ((1.0/12.0))*H_B0(0,-2,0)*u1_B0(0,-2,0) + ((2.0/3.0))*H_B0(0,1,0)*u1_B0(0,1,0))*invDelta1block0;

    d1_p_dy = (-(2.0/3.0)*p_B0(0,-1,0) - (1.0/12.0)*p_B0(0,2,0) + ((1.0/12.0))*p_B0(0,-2,0) +
      ((2.0/3.0))*p_B0(0,1,0))*invDelta1block0;

    d1_rho_dy = (-(2.0/3.0)*rho_B0(0,-1,0) - (1.0/12.0)*rho_B0(0,2,0) + ((1.0/12.0))*rho_B0(0,-2,0) +
      ((2.0/3.0))*rho_B0(0,1,0))*invDelta1block0;

    d1_rhou0_dy = (-(2.0/3.0)*rhou0_B0(0,-1,0) - (1.0/12.0)*rhou0_B0(0,2,0) + ((1.0/12.0))*rhou0_B0(0,-2,0) +
      ((2.0/3.0))*rhou0_B0(0,1,0))*invDelta1block0;

    d1_rhou1_dy = (-(2.0/3.0)*rhou1_B0(0,-1,0) - (1.0/12.0)*rhou1_B0(0,2,0) + ((1.0/12.0))*rhou1_B0(0,-2,0) +
      ((2.0/3.0))*rhou1_B0(0,1,0))*invDelta1block0;

    d1_rhou1u0_dy = (-(2.0/3.0)*u0_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0)*rhou1_B0(0,2,0) +
      ((1.0/12.0))*u0_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u0_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

    d1_rhou1u1_dy = (-(2.0/3.0)*u1_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0)*rhou1_B0(0,2,0) +
      ((1.0/12.0))*u1_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u1_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

    d1_rhou1u2_dy = (-(2.0/3.0)*u2_B0(0,-1,0)*rhou1_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0)*rhou1_B0(0,2,0) +
      ((1.0/12.0))*u2_B0(0,-2,0)*rhou1_B0(0,-2,0) + ((2.0/3.0))*u2_B0(0,1,0)*rhou1_B0(0,1,0))*invDelta1block0;

    d1_rhou2_dy = (-(2.0/3.0)*rhou2_B0(0,-1,0) - (1.0/12.0)*rhou2_B0(0,2,0) + ((1.0/12.0))*rhou2_B0(0,-2,0) +
      ((2.0/3.0))*rhou2_B0(0,1,0))*invDelta1block0;

    d1_u0u1_dy = (-(2.0/3.0)*u0_B0(0,-1,0)*u1_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0)*u1_B0(0,2,0) +
      ((1.0/12.0))*u0_B0(0,-2,0)*u1_B0(0,-2,0) + ((2.0/3.0))*u0_B0(0,1,0)*u1_B0(0,1,0))*invDelta1block0;

    d1_u1u1_dy = (-(2.0/3.0)*(u1_B0(0,-1,0)*u1_B0(0,-1,0)) - (1.0/12.0)*(u1_B0(0,2,0)*u1_B0(0,2,0)) +
      ((1.0/12.0))*(u1_B0(0,-2,0)*u1_B0(0,-2,0)) + ((2.0/3.0))*(u1_B0(0,1,0)*u1_B0(0,1,0)))*invDelta1block0;

    d1_u1u2_dy = (-(2.0/3.0)*u1_B0(0,-1,0)*u2_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0)*u2_B0(0,2,0) +
      ((1.0/12.0))*u1_B0(0,-2,0)*u2_B0(0,-2,0) + ((2.0/3.0))*u1_B0(0,1,0)*u2_B0(0,1,0))*invDelta1block0;

    d1_H_dz = (-(2.0/3.0)*H_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2) + ((1.0/12.0))*H_B0(0,0,-2) +
      ((2.0/3.0))*H_B0(0,0,1))*invDelta2block0;

    d1_Hrho_dz = (-(2.0/3.0)*H_B0(0,0,-1)*rho_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*rho_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*rho_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*rho_B0(0,0,1))*invDelta2block0;

    d1_Hrhou2_dz = (-(2.0/3.0)*H_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_Hu2_dz = (-(2.0/3.0)*H_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*H_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*H_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*H_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_p_dz = (-(2.0/3.0)*p_B0(0,0,-1) - (1.0/12.0)*p_B0(0,0,2) + ((1.0/12.0))*p_B0(0,0,-2) +
      ((2.0/3.0))*p_B0(0,0,1))*invDelta2block0;

    d1_rho_dz = (-(2.0/3.0)*rho_B0(0,0,-1) - (1.0/12.0)*rho_B0(0,0,2) + ((1.0/12.0))*rho_B0(0,0,-2) +
      ((2.0/3.0))*rho_B0(0,0,1))*invDelta2block0;

    d1_rhou0_dz = (-(2.0/3.0)*rhou0_B0(0,0,-1) - (1.0/12.0)*rhou0_B0(0,0,2) + ((1.0/12.0))*rhou0_B0(0,0,-2) +
      ((2.0/3.0))*rhou0_B0(0,0,1))*invDelta2block0;

    d1_rhou1_dz = (-(2.0/3.0)*rhou1_B0(0,0,-1) - (1.0/12.0)*rhou1_B0(0,0,2) + ((1.0/12.0))*rhou1_B0(0,0,-2) +
      ((2.0/3.0))*rhou1_B0(0,0,1))*invDelta2block0;

    d1_rhou2_dz = (-(2.0/3.0)*rhou2_B0(0,0,-1) - (1.0/12.0)*rhou2_B0(0,0,2) + ((1.0/12.0))*rhou2_B0(0,0,-2) +
      ((2.0/3.0))*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u0_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u0_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u1_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u1_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_rhou2u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1)*rhou2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2)*rhou2_B0(0,0,2) +
      ((1.0/12.0))*u2_B0(0,0,-2)*rhou2_B0(0,0,-2) + ((2.0/3.0))*u2_B0(0,0,1)*rhou2_B0(0,0,1))*invDelta2block0;

    d1_u0u2_dz = (-(2.0/3.0)*u0_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*u0_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*u0_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_u1u2_dz = (-(2.0/3.0)*u1_B0(0,0,-1)*u2_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2)*u2_B0(0,0,2) +
      ((1.0/12.0))*u1_B0(0,0,-2)*u2_B0(0,0,-2) + ((2.0/3.0))*u1_B0(0,0,1)*u2_B0(0,0,1))*invDelta2block0;

    d1_u2u2_dz = (-(2.0/3.0)*(u2_B0(0,0,-1)*u2_B0(0,0,-1)) - (1.0/12.0)*(u2_B0(0,0,2)*u2_B0(0,0,2)) +
      ((1.0/12.0))*(u2_B0(0,0,-2)*u2_B0(0,0,-2)) + ((2.0/3.0))*(u2_B0(0,0,1)*u2_B0(0,0,1)))*invDelta2block0;

    Residual0_B0(0,0,0) = -(1.0/2.0)*d1_rhou0_dx - (1.0/2.0)*d1_rhou1_dy - (1.0/2.0)*d1_rhou2_dz -
      (1.0/2.0)*(wk0_B0(0,0,0) + wk4_B0(0,0,0) + wk8_B0(0,0,0))*rho_B0(0,0,0) - (1.0/2.0)*u0_B0(0,0,0)*d1_rho_dx -
      (1.0/2.0)*u1_B0(0,0,0)*d1_rho_dy - (1.0/2.0)*u2_B0(0,0,0)*d1_rho_dz;

    Residual1_B0(0,0,0) = -d1_p_dx - (1.0/4.0)*d1_rhou0u0_dx - (1.0/4.0)*d1_rhou1u0_dy - (1.0/4.0)*d1_rhou2u0_dz -
      (1.0/2.0)*u0_B0(0,0,0)*d1_rhou0_dx - (1.0/2.0)*wk0_B0(0,0,0)*rhou0_B0(0,0,0) -
      (1.0/4.0)*(u0_B0(0,0,0)*u0_B0(0,0,0))*d1_rho_dx - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou1_dy -
      (1.0/4.0)*u0_B0(0,0,0)*d1_rhou2_dz - (1.0/4.0)*u1_B0(0,0,0)*d1_rhou0_dy - (1.0/4.0)*u2_B0(0,0,0)*d1_rhou0_dz -
      (1.0/4.0)*rho_B0(0,0,0)*d1_u0u0_dx - (1.0/4.0)*rho_B0(0,0,0)*d1_u0u1_dy - (1.0/4.0)*rho_B0(0,0,0)*d1_u0u2_dz -
      (1.0/4.0)*wk3_B0(0,0,0)*rhou1_B0(0,0,0) - (1.0/4.0)*wk4_B0(0,0,0)*rhou0_B0(0,0,0) -
      (1.0/4.0)*wk6_B0(0,0,0)*rhou2_B0(0,0,0) - (1.0/4.0)*wk8_B0(0,0,0)*rhou0_B0(0,0,0) -
      (1.0/4.0)*u0_B0(0,0,0)*u1_B0(0,0,0)*d1_rho_dy - (1.0/4.0)*u0_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz;

    Residual2_B0(0,0,0) = -d1_p_dy - (1.0/4.0)*d1_rhou0u1_dx - (1.0/4.0)*d1_rhou1u1_dy - (1.0/4.0)*d1_rhou2u1_dz -
      (1.0/4.0)*(u1_B0(0,0,0)*u1_B0(0,0,0))*d1_rho_dy - (1.0/4.0)*(wk0_B0(0,0,0) + wk4_B0(0,0,0) +
      wk8_B0(0,0,0))*rhou1_B0(0,0,0) - (1.0/4.0)*(d1_rhou0_dx + d1_rhou1_dy + d1_rhou2_dz)*u1_B0(0,0,0) -
      (1.0/4.0)*(d1_u0u1_dx + d1_u1u1_dy + d1_u1u2_dz)*rho_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou1_dx -
      (1.0/4.0)*u1_B0(0,0,0)*d1_rhou1_dy - (1.0/4.0)*u2_B0(0,0,0)*d1_rhou1_dz - (1.0/4.0)*wk1_B0(0,0,0)*rhou0_B0(0,0,0)
      - (1.0/4.0)*wk4_B0(0,0,0)*rhou1_B0(0,0,0) - (1.0/4.0)*wk7_B0(0,0,0)*rhou2_B0(0,0,0) -
      (1.0/4.0)*u0_B0(0,0,0)*u1_B0(0,0,0)*d1_rho_dx - (1.0/4.0)*u1_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz;

    Residual3_B0(0,0,0) = -d1_p_dz - (1.0/4.0)*d1_rhou0u2_dx - (1.0/4.0)*d1_rhou1u2_dy - (1.0/4.0)*d1_rhou2u2_dz -
      (1.0/4.0)*(u2_B0(0,0,0)*u2_B0(0,0,0))*d1_rho_dz - (1.0/4.0)*(wk0_B0(0,0,0) + wk4_B0(0,0,0) +
      wk8_B0(0,0,0))*rhou2_B0(0,0,0) - (1.0/4.0)*(d1_rhou0_dx + d1_rhou1_dy + d1_rhou2_dz)*u2_B0(0,0,0) -
      (1.0/4.0)*(d1_u0u2_dx + d1_u1u2_dy + d1_u2u2_dz)*rho_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_rhou2_dx -
      (1.0/4.0)*u1_B0(0,0,0)*d1_rhou2_dy - (1.0/4.0)*u2_B0(0,0,0)*d1_rhou2_dz - (1.0/4.0)*wk2_B0(0,0,0)*rhou0_B0(0,0,0)
      - (1.0/4.0)*wk5_B0(0,0,0)*rhou1_B0(0,0,0) - (1.0/4.0)*wk8_B0(0,0,0)*rhou2_B0(0,0,0) -
      (1.0/4.0)*u0_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dx - (1.0/4.0)*u1_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dy;

    Residual4_B0(0,0,0) = -(1.0/4.0)*d1_Hrhou0_dx - (1.0/4.0)*d1_Hrhou1_dy - (1.0/4.0)*d1_Hrhou2_dz -
      (1.0/4.0)*(d1_Hu0_dx + d1_Hu1_dy + d1_Hu2_dz)*rho_B0(0,0,0) - (1.0/4.0)*(d1_rhou0_dx + d1_rhou1_dy +
      d1_rhou2_dz)*H_B0(0,0,0) - (1.0/4.0)*u0_B0(0,0,0)*d1_Hrho_dx - (1.0/4.0)*u1_B0(0,0,0)*d1_Hrho_dy -
      (1.0/4.0)*u2_B0(0,0,0)*d1_Hrho_dz - (1.0/4.0)*rhou0_B0(0,0,0)*d1_H_dx - (1.0/4.0)*rhou1_B0(0,0,0)*d1_H_dy -
      (1.0/4.0)*rhou2_B0(0,0,0)*d1_H_dz - (1.0/4.0)*(wk0_B0(0,0,0) + wk4_B0(0,0,0) +
      wk8_B0(0,0,0))*H_B0(0,0,0)*rho_B0(0,0,0) - (1.0/4.0)*H_B0(0,0,0)*u0_B0(0,0,0)*d1_rho_dx -
      (1.0/4.0)*H_B0(0,0,0)*u1_B0(0,0,0)*d1_rho_dy - (1.0/4.0)*H_B0(0,0,0)*u2_B0(0,0,0)*d1_rho_dz;

}

 void opensbliblock00Kernel030(const ACC<double> &T_B0, const ACC<double> &mu_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, const ACC<double> &wk0_B0, const ACC<double> &wk1_B0, const ACC<double>
&wk2_B0, const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, const ACC<double> &wk5_B0, const ACC<double> &wk6_B0,
const ACC<double> &wk7_B0, const ACC<double> &wk8_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double>
&Residual3_B0, ACC<double> &Residual4_B0)
{
   double d1_T_dx = 0.0;
   double d1_T_dy = 0.0;
   double d1_T_dz = 0.0;
   double d1_mu_dx = 0.0;
   double d1_mu_dy = 0.0;
   double d1_mu_dz = 0.0;
   double d1_wk0_dy = 0.0;
   double d1_wk0_dz = 0.0;
   double d1_wk1_dy = 0.0;
   double d1_wk2_dz = 0.0;
   double d1_wk4_dz = 0.0;
   double d1_wk5_dz = 0.0;
   double d2_T_dx = 0.0;
   double d2_T_dy = 0.0;
   double d2_T_dz = 0.0;
   double d2_u0_dx = 0.0;
   double d2_u0_dy = 0.0;
   double d2_u0_dz = 0.0;
   double d2_u1_dx = 0.0;
   double d2_u1_dy = 0.0;
   double d2_u1_dz = 0.0;
   double d2_u2_dx = 0.0;
   double d2_u2_dy = 0.0;
   double d2_u2_dz = 0.0;
    d1_T_dx = (-(2.0/3.0)*T_B0(-1,0,0) - (1.0/12.0)*T_B0(2,0,0) + ((1.0/12.0))*T_B0(-2,0,0) +
      ((2.0/3.0))*T_B0(1,0,0))*invDelta0block0;

    d2_T_dx = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(-2,0,0) - (1.0/12.0)*T_B0(2,0,0) + ((4.0/3.0))*T_B0(1,0,0) +
      ((4.0/3.0))*T_B0(-1,0,0))*inv2Delta0block0;

    d1_mu_dx = (-(2.0/3.0)*mu_B0(-1,0,0) - (1.0/12.0)*mu_B0(2,0,0) + ((1.0/12.0))*mu_B0(-2,0,0) +
      ((2.0/3.0))*mu_B0(1,0,0))*invDelta0block0;

    d2_u0_dx = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(-2,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((4.0/3.0))*u0_B0(1,0,0)
      + ((4.0/3.0))*u0_B0(-1,0,0))*inv2Delta0block0;

    d2_u1_dx = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(-2,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((4.0/3.0))*u1_B0(1,0,0)
      + ((4.0/3.0))*u1_B0(-1,0,0))*inv2Delta0block0;

    d2_u2_dx = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(-2,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((4.0/3.0))*u2_B0(1,0,0)
      + ((4.0/3.0))*u2_B0(-1,0,0))*inv2Delta0block0;

    d1_T_dy = (-(2.0/3.0)*T_B0(0,-1,0) - (1.0/12.0)*T_B0(0,2,0) + ((1.0/12.0))*T_B0(0,-2,0) +
      ((2.0/3.0))*T_B0(0,1,0))*invDelta1block0;

    d2_T_dy = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,-2,0) - (1.0/12.0)*T_B0(0,2,0) + ((4.0/3.0))*T_B0(0,1,0) +
      ((4.0/3.0))*T_B0(0,-1,0))*inv2Delta1block0;

    d1_mu_dy = (-(2.0/3.0)*mu_B0(0,-1,0) - (1.0/12.0)*mu_B0(0,2,0) + ((1.0/12.0))*mu_B0(0,-2,0) +
      ((2.0/3.0))*mu_B0(0,1,0))*invDelta1block0;

    d2_u0_dy = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,-2,0) - (1.0/12.0)*u0_B0(0,2,0) + ((4.0/3.0))*u0_B0(0,1,0)
      + ((4.0/3.0))*u0_B0(0,-1,0))*inv2Delta1block0;

    d2_u1_dy = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,-2,0) - (1.0/12.0)*u1_B0(0,2,0) + ((4.0/3.0))*u1_B0(0,1,0)
      + ((4.0/3.0))*u1_B0(0,-1,0))*inv2Delta1block0;

    d2_u2_dy = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,-2,0) - (1.0/12.0)*u2_B0(0,2,0) + ((4.0/3.0))*u2_B0(0,1,0)
      + ((4.0/3.0))*u2_B0(0,-1,0))*inv2Delta1block0;

    d1_wk0_dy = (-(2.0/3.0)*wk0_B0(0,-1,0) - (1.0/12.0)*wk0_B0(0,2,0) + ((1.0/12.0))*wk0_B0(0,-2,0) +
      ((2.0/3.0))*wk0_B0(0,1,0))*invDelta1block0;

    d1_wk1_dy = (-(2.0/3.0)*wk1_B0(0,-1,0) - (1.0/12.0)*wk1_B0(0,2,0) + ((1.0/12.0))*wk1_B0(0,-2,0) +
      ((2.0/3.0))*wk1_B0(0,1,0))*invDelta1block0;

    d1_T_dz = (-(2.0/3.0)*T_B0(0,0,-1) - (1.0/12.0)*T_B0(0,0,2) + ((1.0/12.0))*T_B0(0,0,-2) +
      ((2.0/3.0))*T_B0(0,0,1))*invDelta2block0;

    d2_T_dz = (-(5.0/2.0)*T_B0(0,0,0) - (1.0/12.0)*T_B0(0,0,-2) - (1.0/12.0)*T_B0(0,0,2) + ((4.0/3.0))*T_B0(0,0,1) +
      ((4.0/3.0))*T_B0(0,0,-1))*inv2Delta2block0;

    d1_mu_dz = (-(2.0/3.0)*mu_B0(0,0,-1) - (1.0/12.0)*mu_B0(0,0,2) + ((1.0/12.0))*mu_B0(0,0,-2) +
      ((2.0/3.0))*mu_B0(0,0,1))*invDelta2block0;

    d2_u0_dz = (-(5.0/2.0)*u0_B0(0,0,0) - (1.0/12.0)*u0_B0(0,0,-2) - (1.0/12.0)*u0_B0(0,0,2) + ((4.0/3.0))*u0_B0(0,0,1)
      + ((4.0/3.0))*u0_B0(0,0,-1))*inv2Delta2block0;

    d2_u1_dz = (-(5.0/2.0)*u1_B0(0,0,0) - (1.0/12.0)*u1_B0(0,0,-2) - (1.0/12.0)*u1_B0(0,0,2) + ((4.0/3.0))*u1_B0(0,0,1)
      + ((4.0/3.0))*u1_B0(0,0,-1))*inv2Delta2block0;

    d2_u2_dz = (-(5.0/2.0)*u2_B0(0,0,0) - (1.0/12.0)*u2_B0(0,0,-2) - (1.0/12.0)*u2_B0(0,0,2) + ((4.0/3.0))*u2_B0(0,0,1)
      + ((4.0/3.0))*u2_B0(0,0,-1))*inv2Delta2block0;

    d1_wk0_dz = (-(2.0/3.0)*wk0_B0(0,0,-1) - (1.0/12.0)*wk0_B0(0,0,2) + ((1.0/12.0))*wk0_B0(0,0,-2) +
      ((2.0/3.0))*wk0_B0(0,0,1))*invDelta2block0;

    d1_wk2_dz = (-(2.0/3.0)*wk2_B0(0,0,-1) - (1.0/12.0)*wk2_B0(0,0,2) + ((1.0/12.0))*wk2_B0(0,0,-2) +
      ((2.0/3.0))*wk2_B0(0,0,1))*invDelta2block0;

    d1_wk4_dz = (-(2.0/3.0)*wk4_B0(0,0,-1) - (1.0/12.0)*wk4_B0(0,0,2) + ((1.0/12.0))*wk4_B0(0,0,-2) +
      ((2.0/3.0))*wk4_B0(0,0,1))*invDelta2block0;

    d1_wk5_dz = (-(2.0/3.0)*wk5_B0(0,0,-1) - (1.0/12.0)*wk5_B0(0,0,2) + ((1.0/12.0))*wk5_B0(0,0,-2) +
      ((2.0/3.0))*wk5_B0(0,0,1))*invDelta2block0;

    Residual1_B0(0,0,0) = (wk1_B0(0,0,0) + wk3_B0(0,0,0))*invRe*d1_mu_dy + (wk2_B0(0,0,0) +
      wk6_B0(0,0,0))*invRe*d1_mu_dz + (-(2.0/3.0)*wk4_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) +
      ((4.0/3.0))*wk0_B0(0,0,0))*invRe*d1_mu_dx + (((1.0/3.0))*d1_wk1_dy + ((1.0/3.0))*d1_wk2_dz + ((4.0/3.0))*d2_u0_dx
      + d2_u0_dy + d2_u0_dz)*invRe*mu_B0(0,0,0) + Residual1_B0(0,0,0);

    Residual2_B0(0,0,0) = (wk1_B0(0,0,0) + wk3_B0(0,0,0))*invRe*d1_mu_dx + (wk5_B0(0,0,0) +
      wk7_B0(0,0,0))*invRe*d1_mu_dz + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) +
      ((4.0/3.0))*wk4_B0(0,0,0))*invRe*d1_mu_dy + (((1.0/3.0))*d1_wk0_dy + ((1.0/3.0))*d1_wk5_dz + ((4.0/3.0))*d2_u1_dy
      + d2_u1_dx + d2_u1_dz)*invRe*mu_B0(0,0,0) + Residual2_B0(0,0,0);

    Residual3_B0(0,0,0) = (wk2_B0(0,0,0) + wk6_B0(0,0,0))*invRe*d1_mu_dx + (wk5_B0(0,0,0) +
      wk7_B0(0,0,0))*invRe*d1_mu_dy + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk4_B0(0,0,0) +
      ((4.0/3.0))*wk8_B0(0,0,0))*invRe*d1_mu_dz + (((1.0/3.0))*d1_wk0_dz + ((1.0/3.0))*d1_wk4_dz + ((4.0/3.0))*d2_u2_dz
      + d2_u2_dx + d2_u2_dy)*invRe*mu_B0(0,0,0) + Residual3_B0(0,0,0);

    Residual4_B0(0,0,0) = (wk1_B0(0,0,0) + wk3_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk1_B0(0,0,0) + (wk1_B0(0,0,0) +
      wk3_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk3_B0(0,0,0) + (wk1_B0(0,0,0) + wk3_B0(0,0,0))*invRe*u0_B0(0,0,0)*d1_mu_dy +
      (wk1_B0(0,0,0) + wk3_B0(0,0,0))*invRe*u1_B0(0,0,0)*d1_mu_dx + (wk2_B0(0,0,0) +
      wk6_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk2_B0(0,0,0) + (wk2_B0(0,0,0) + wk6_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk6_B0(0,0,0)
      + (wk2_B0(0,0,0) + wk6_B0(0,0,0))*invRe*u0_B0(0,0,0)*d1_mu_dz + (wk2_B0(0,0,0) +
      wk6_B0(0,0,0))*invRe*u2_B0(0,0,0)*d1_mu_dx + (wk5_B0(0,0,0) + wk7_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk5_B0(0,0,0) +
      (wk5_B0(0,0,0) + wk7_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk7_B0(0,0,0) + (wk5_B0(0,0,0) +
      wk7_B0(0,0,0))*invRe*u1_B0(0,0,0)*d1_mu_dz + (wk5_B0(0,0,0) + wk7_B0(0,0,0))*invRe*u2_B0(0,0,0)*d1_mu_dy +
      (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk4_B0(0,0,0) + ((4.0/3.0))*wk8_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk8_B0(0,0,0)
      + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk4_B0(0,0,0) + ((4.0/3.0))*wk8_B0(0,0,0))*invRe*u2_B0(0,0,0)*d1_mu_dz +
      (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) + ((4.0/3.0))*wk4_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk4_B0(0,0,0)
      + (-(2.0/3.0)*wk0_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) + ((4.0/3.0))*wk4_B0(0,0,0))*invRe*u1_B0(0,0,0)*d1_mu_dy +
      (-(2.0/3.0)*wk4_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) + ((4.0/3.0))*wk0_B0(0,0,0))*invRe*mu_B0(0,0,0)*wk0_B0(0,0,0)
      + (-(2.0/3.0)*wk4_B0(0,0,0) - (2.0/3.0)*wk8_B0(0,0,0) + ((4.0/3.0))*wk0_B0(0,0,0))*invRe*u0_B0(0,0,0)*d1_mu_dx +
      (((1.0/3.0))*d1_wk0_dy + ((1.0/3.0))*d1_wk5_dz + ((4.0/3.0))*d2_u1_dy + d2_u1_dx +
      d2_u1_dz)*invRe*mu_B0(0,0,0)*u1_B0(0,0,0) + (((1.0/3.0))*d1_wk0_dz + ((1.0/3.0))*d1_wk4_dz + ((4.0/3.0))*d2_u2_dz
      + d2_u2_dx + d2_u2_dy)*invRe*mu_B0(0,0,0)*u2_B0(0,0,0) + (((1.0/3.0))*d1_wk1_dy + ((1.0/3.0))*d1_wk2_dz +
      ((4.0/3.0))*d2_u0_dx + d2_u0_dy + d2_u0_dz)*invRe*mu_B0(0,0,0)*u0_B0(0,0,0) + (d2_T_dx + d2_T_dy +
      d2_T_dz)*invPr*invRe*inv2Minf*inv_gamma_m1*mu_B0(0,0,0) + invPr*invRe*inv2Minf*inv_gamma_m1*d1_T_dx*d1_mu_dx +
      invPr*invRe*inv2Minf*inv_gamma_m1*d1_T_dy*d1_mu_dy + invPr*invRe*inv2Minf*inv_gamma_m1*d1_T_dz*d1_mu_dz +
      Residual4_B0(0,0,0);

}

 void opensbliblock00Kernel050(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
&Residual2_B0, const ACC<double> &Residual3_B0, const ACC<double> &Residual4_B0, ACC<double> &rhoE_B0, ACC<double>
&rhoE_RKold_B0, ACC<double> &rho_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_B0, ACC<double> &rhou0_RKold_B0,
ACC<double> &rhou1_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_B0, ACC<double> &rhou2_RKold_B0, const double
*rkA, const double *rkB)
{
   rho_RKold_B0(0,0,0) = rkA[0]*rho_RKold_B0(0,0,0) + dt*Residual0_B0(0,0,0);

   rho_B0(0,0,0) = rkB[0]*rho_RKold_B0(0,0,0) + rho_B0(0,0,0);

   rhou0_RKold_B0(0,0,0) = rkA[0]*rhou0_RKold_B0(0,0,0) + dt*Residual1_B0(0,0,0);

   rhou0_B0(0,0,0) = rkB[0]*rhou0_RKold_B0(0,0,0) + rhou0_B0(0,0,0);

   rhou1_RKold_B0(0,0,0) = rkA[0]*rhou1_RKold_B0(0,0,0) + dt*Residual2_B0(0,0,0);

   rhou1_B0(0,0,0) = rkB[0]*rhou1_RKold_B0(0,0,0) + rhou1_B0(0,0,0);

   rhou2_RKold_B0(0,0,0) = rkA[0]*rhou2_RKold_B0(0,0,0) + dt*Residual3_B0(0,0,0);

   rhou2_B0(0,0,0) = rkB[0]*rhou2_RKold_B0(0,0,0) + rhou2_B0(0,0,0);

   rhoE_RKold_B0(0,0,0) = rkA[0]*rhoE_RKold_B0(0,0,0) + dt*Residual4_B0(0,0,0);

   rhoE_B0(0,0,0) = rkB[0]*rhoE_RKold_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel043(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const
ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, ACC<double> &a_B0, ACC<double> &u0_B0, ACC<double> &u1_B0,
ACC<double> &u2_B0, ACC<double> &p_B0)
{
   double inv_rho = 0.0;
   inv_rho = 1.0/rho_B0(0,0,0);

   u0_B0(0,0,0) = rhou0_B0(0,0,0)*inv_rho;

   u1_B0(0,0,0) = rhou1_B0(0,0,0)*inv_rho;

   u2_B0(0,0,0) = rhou2_B0(0,0,0)*inv_rho;

    p_B0(0,0,0) = (-1 + gama)*(-(0.5*(rhou0_B0(0,0,0)*rhou0_B0(0,0,0)) + 0.5*(rhou1_B0(0,0,0)*rhou1_B0(0,0,0)) +
      0.5*(rhou2_B0(0,0,0)*rhou2_B0(0,0,0)))*inv_rho + rhoE_B0(0,0,0));

   a_B0(0,0,0) = sqrt((gama*p_B0(0,0,0)*inv_rho));

}

 void opensbliblock00Kernel044(const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&kappa_B0)
{
   double d1_u0_dx = 0.0;
   double d1_u0_dy = 0.0;
   double d1_u0_dz = 0.0;
   double d1_u1_dx = 0.0;
   double d1_u1_dy = 0.0;
   double d1_u1_dz = 0.0;
   double d1_u2_dx = 0.0;
   double d1_u2_dy = 0.0;
   double d1_u2_dz = 0.0;
    d1_u0_dy = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
      ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

    d1_u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

    d1_u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

    d1_u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

    d1_u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

    d1_u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

    d1_u2_dy = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
      ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

    d1_u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

    d1_u1_dy = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
      ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

    kappa_B0(0,0,0) = ((d1_u0_dx + d1_u1_dy + d1_u2_dz)*(d1_u0_dx + d1_u1_dy + d1_u2_dz))*(0.5 - 0.5*tanh(2.5 +
      500.0*(d1_u0_dx + d1_u1_dy + d1_u2_dz)/sqrt(((Delta0block0*Delta0block0) + (Delta1block0*Delta1block0) +
      (Delta2block0*Delta2block0)))))/(1.0e-40 + ((-d1_u0_dy + d1_u1_dx)*(-d1_u0_dy + d1_u1_dx)) + ((-d1_u1_dz +
      d1_u2_dy)*(-d1_u1_dz + d1_u2_dy)) + ((-d1_u2_dx + d1_u0_dz)*(-d1_u2_dx + d1_u0_dz)) + ((d1_u0_dx + d1_u1_dy +
      d1_u2_dz)*(d1_u0_dx + d1_u1_dy + d1_u2_dz)));

}

 void opensbliblock00Kernel045(ACC<double> &Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0,
ACC<double> &Residual3_B0, ACC<double> &Residual4_B0, ACC<double> &rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double>
&rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double> &rhou2_RKold_B0, ACC<double> &wk0_B0, ACC<double> &wk1_B0,
ACC<double> &wk2_B0, ACC<double> &wk3_B0, ACC<double> &wk4_B0)
{
   wk0_B0(0,0,0) = 0.0;

   wk1_B0(0,0,0) = 0.0;

   wk2_B0(0,0,0) = 0.0;

   wk3_B0(0,0,0) = 0.0;

   wk4_B0(0,0,0) = 0.0;

   Residual0_B0(0,0,0) = 0.0;

   Residual1_B0(0,0,0) = 0.0;

   Residual2_B0(0,0,0) = 0.0;

   Residual3_B0(0,0,0) = 0.0;

   Residual4_B0(0,0,0) = 0.0;

   rho_RKold_B0(0,0,0) = 0.0;

   rhou0_RKold_B0(0,0,0) = 0.0;

   rhou1_RKold_B0(0,0,0) = 0.0;

   rhou2_RKold_B0(0,0,0) = 0.0;

   rhoE_RKold_B0(0,0,0) = 0.0;

}

 void opensbliblock00Kernel046(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&wk0_B0, ACC<double> &wk1_B0, ACC<double> &wk2_B0, ACC<double> &wk3_B0, ACC<double> &wk4_B0)
{
   double AVG_0_0_LEV_00 = 0.0;
   double AVG_0_0_LEV_01 = 0.0;
   double AVG_0_0_LEV_02 = 0.0;
   double AVG_0_0_LEV_03 = 0.0;
   double AVG_0_0_LEV_04 = 0.0;
   double AVG_0_0_LEV_10 = 0.0;
   double AVG_0_0_LEV_13 = 0.0;
   double AVG_0_0_LEV_20 = 0.0;
   double AVG_0_0_LEV_22 = 0.0;
   double AVG_0_0_LEV_30 = 0.0;
   double AVG_0_0_LEV_31 = 0.0;
   double AVG_0_0_LEV_32 = 0.0;
   double AVG_0_0_LEV_33 = 0.0;
   double AVG_0_0_LEV_34 = 0.0;
   double AVG_0_0_LEV_40 = 0.0;
   double AVG_0_0_LEV_41 = 0.0;
   double AVG_0_0_LEV_42 = 0.0;
   double AVG_0_0_LEV_43 = 0.0;
   double AVG_0_0_LEV_44 = 0.0;
   double AVG_0_a = 0.0;
   double AVG_0_inv_rho = 0.0;
   double AVG_0_rho = 0.0;
   double AVG_0_u0 = 0.0;
   double AVG_0_u1 = 0.0;
   double AVG_0_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(-2,0,0), fmax(kappa_B0(0,0,0), fmax(kappa_B0(-1,0,0), fmax(kappa_B0(2,0,0), fmax(kappa_B0(-3,0,0),
      kappa_B0(1,0,0)))))) > Ducros_check){

      AVG_0_rho = sqrt((rho_B0(0,0,0)*rho_B0(1,0,0)));

      AVG_0_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(1,0,0))));

      AVG_0_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u0_B0(1,0,0))*AVG_0_inv_rho;

      AVG_0_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u1_B0(1,0,0))*AVG_0_inv_rho;

      AVG_0_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(1,0,0))*u2_B0(1,0,0))*AVG_0_inv_rho;

       AVG_0_a = sqrt(((-(1.0/2.0)*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) + (AVG_0_u2*AVG_0_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(1,0,0) +
            rhoE_B0(1,0,0))/sqrt(rho_B0(1,0,0)))*AVG_0_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_0_a);

      inv_AVG_rho = 1.0/(AVG_0_rho);

       AVG_0_0_LEV_00 = -(1.0/2.0)*(-2 - (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_0_u2*AVG_0_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_0_u0*AVG_0_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_0_u1*AVG_0_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_0_u2*AVG_0_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_0_0_LEV_01 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u0;

      AVG_0_0_LEV_02 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u1;

      AVG_0_0_LEV_03 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_0_u2;

      AVG_0_0_LEV_04 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

      AVG_0_0_LEV_10 = -AVG_0_u2*inv_AVG_rho;

      AVG_0_0_LEV_13 = inv_AVG_rho;

      AVG_0_0_LEV_20 = AVG_0_u1*inv_AVG_rho;

      AVG_0_0_LEV_22 = -inv_AVG_rho;

       AVG_0_0_LEV_30 = -0.353553390593274*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) + (AVG_0_u2*AVG_0_u2) -
            (AVG_0_u0*AVG_0_u0)*gama - (AVG_0_u1*AVG_0_u1)*gama - (AVG_0_u2*AVG_0_u2)*gama +
            2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_31 = 0.707106781186547*(-gama*AVG_0_u0 + AVG_0_a + AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_32 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_33 = -0.707106781186547*gamma_m1*AVG_0_u2*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_0_0_LEV_40 = 0.353553390593274*(-(AVG_0_u0*AVG_0_u0) - (AVG_0_u1*AVG_0_u1) - (AVG_0_u2*AVG_0_u2) +
            (AVG_0_u0*AVG_0_u0)*gama + (AVG_0_u1*AVG_0_u1)*gama + (AVG_0_u2*AVG_0_u2)*gama +
            2*AVG_0_a*AVG_0_u0)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_41 = -0.707106781186547*(-AVG_0_u0 + gama*AVG_0_u0 + AVG_0_a)*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_42 = -0.707106781186547*gamma_m1*AVG_0_u1*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_43 = -0.707106781186547*gamma_m1*AVG_0_u2*inv_AVG_a*inv_AVG_rho;

      AVG_0_0_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       CF_00 = p_B0(-1,0,0)*AVG_0_0_LEV_01 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_00 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_04
            + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 +
            u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

      CF_10 = rhou0_B0(-1,0,0)*AVG_0_0_LEV_10 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_13;

      CF_20 = rhou0_B0(-1,0,0)*AVG_0_0_LEV_20 + u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_22;

       CF_30 = p_B0(-1,0,0)*AVG_0_0_LEV_31 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_30 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_34
            + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 +
            u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

       CF_40 = p_B0(-1,0,0)*AVG_0_0_LEV_41 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_40 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_44
            + u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 +
            u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

       CS_00 = rho_B0(-1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

      CS_10 = rho_B0(-1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_13;

      CS_20 = rho_B0(-1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-1,0,0)*AVG_0_0_LEV_22;

       CS_30 = rho_B0(-1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

       CS_40 = rho_B0(-1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

       CF_01 = p_B0(0,0,0)*AVG_0_0_LEV_01 + rhou0_B0(0,0,0)*AVG_0_0_LEV_00 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_04 +
            u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_04 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_01 +
            u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_02 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

      CF_11 = rhou0_B0(0,0,0)*AVG_0_0_LEV_10 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_13;

      CF_21 = rhou0_B0(0,0,0)*AVG_0_0_LEV_20 + u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_22;

       CF_31 = p_B0(0,0,0)*AVG_0_0_LEV_31 + rhou0_B0(0,0,0)*AVG_0_0_LEV_30 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_34 +
            u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_34 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_31 +
            u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_32 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

       CF_41 = p_B0(0,0,0)*AVG_0_0_LEV_41 + rhou0_B0(0,0,0)*AVG_0_0_LEV_40 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_44 +
            u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_44 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_41 +
            u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_42 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

       CS_01 = rho_B0(0,0,0)*AVG_0_0_LEV_00 + rhoE_B0(0,0,0)*AVG_0_0_LEV_04 + rhou0_B0(0,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_02 + rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

      CS_11 = rho_B0(0,0,0)*AVG_0_0_LEV_10 + rhou2_B0(0,0,0)*AVG_0_0_LEV_13;

      CS_21 = rho_B0(0,0,0)*AVG_0_0_LEV_20 + rhou1_B0(0,0,0)*AVG_0_0_LEV_22;

       CS_31 = rho_B0(0,0,0)*AVG_0_0_LEV_30 + rhoE_B0(0,0,0)*AVG_0_0_LEV_34 + rhou0_B0(0,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_32 + rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

       CS_41 = rho_B0(0,0,0)*AVG_0_0_LEV_40 + rhoE_B0(0,0,0)*AVG_0_0_LEV_44 + rhou0_B0(0,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(0,0,0)*AVG_0_0_LEV_42 + rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

       CF_02 = p_B0(1,0,0)*AVG_0_0_LEV_01 + rhou0_B0(1,0,0)*AVG_0_0_LEV_00 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_04 +
            u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_04 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_01 +
            u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_02 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

      CF_12 = rhou0_B0(1,0,0)*AVG_0_0_LEV_10 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_13;

      CF_22 = rhou0_B0(1,0,0)*AVG_0_0_LEV_20 + u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_22;

       CF_32 = p_B0(1,0,0)*AVG_0_0_LEV_31 + rhou0_B0(1,0,0)*AVG_0_0_LEV_30 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_34 +
            u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_34 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_31 +
            u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_32 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

       CF_42 = p_B0(1,0,0)*AVG_0_0_LEV_41 + rhou0_B0(1,0,0)*AVG_0_0_LEV_40 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_44 +
            u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_44 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_41 +
            u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_42 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

       CS_02 = rho_B0(1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(1,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

      CS_12 = rho_B0(1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(1,0,0)*AVG_0_0_LEV_13;

      CS_22 = rho_B0(1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(1,0,0)*AVG_0_0_LEV_22;

       CS_32 = rho_B0(1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(1,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

       CS_42 = rho_B0(1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(1,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

       CF_03 = p_B0(2,0,0)*AVG_0_0_LEV_01 + rhou0_B0(2,0,0)*AVG_0_0_LEV_00 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_04 +
            u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_04 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_01 +
            u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_02 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

      CF_13 = rhou0_B0(2,0,0)*AVG_0_0_LEV_10 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_13;

      CF_23 = rhou0_B0(2,0,0)*AVG_0_0_LEV_20 + u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_22;

       CF_33 = p_B0(2,0,0)*AVG_0_0_LEV_31 + rhou0_B0(2,0,0)*AVG_0_0_LEV_30 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_34 +
            u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_34 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_31 +
            u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_32 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

       CF_43 = p_B0(2,0,0)*AVG_0_0_LEV_41 + rhou0_B0(2,0,0)*AVG_0_0_LEV_40 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_44 +
            u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_44 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_41 +
            u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_42 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

       CS_03 = rho_B0(2,0,0)*AVG_0_0_LEV_00 + rhoE_B0(2,0,0)*AVG_0_0_LEV_04 + rhou0_B0(2,0,0)*AVG_0_0_LEV_01 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_02 + rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

      CS_13 = rho_B0(2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(2,0,0)*AVG_0_0_LEV_13;

      CS_23 = rho_B0(2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(2,0,0)*AVG_0_0_LEV_22;

       CS_33 = rho_B0(2,0,0)*AVG_0_0_LEV_30 + rhoE_B0(2,0,0)*AVG_0_0_LEV_34 + rhou0_B0(2,0,0)*AVG_0_0_LEV_31 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_32 + rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

       CS_43 = rho_B0(2,0,0)*AVG_0_0_LEV_40 + rhoE_B0(2,0,0)*AVG_0_0_LEV_44 + rhou0_B0(2,0,0)*AVG_0_0_LEV_41 +
            rhou1_B0(2,0,0)*AVG_0_0_LEV_42 + rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u0_B0(1,0,0)), fabs(u0_B0(0,0,0)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(1,0,0) + u0_B0(1,0,0)), fabs(a_B0(0,0,0) + u0_B0(0,0,0)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u0_B0(0,0,0) + a_B0(0,0,0)), fabs(-u0_B0(1,0,0) + a_B0(1,0,0)));

       beta_0 = ((1.0/4.0))*((CS_02*max_lambda_00 + CF_02)*(CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (CS_02*max_lambda_00 + CF_02))*(CS_01*max_lambda_00
            + CF_01);

       beta_1 = ((1.0/4.0))*((CS_01*max_lambda_00 + CF_01)*(CS_01*max_lambda_00 + CF_01)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - (CS_01*max_lambda_00 + CF_01))*(CS_00*max_lambda_00
            + CF_00);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_0 = (-(1.0/4.0)*(CS_00*max_lambda_00 + CF_00) + ((3.0/4.0))*(CS_01*max_lambda_00 + CF_01))*omega_1 +
            (((1.0/4.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(CS_02*max_lambda_00 + CF_02))*omega_0 + Recon_0;

       beta_0 = ((1.0/4.0))*((-CS_03*max_lambda_00 + CF_03)*(-CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (-CS_03*max_lambda_00 +
            CF_03))*(-CS_02*max_lambda_00 + CF_02);

       beta_1 = ((1.0/4.0))*((-CS_02*max_lambda_00 + CF_02)*(-CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - (-CS_02*max_lambda_00 +
            CF_02))*(-CS_01*max_lambda_00 + CF_01);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = fmax(rj0, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_0 = (-(1.0/4.0)*(-CS_03*max_lambda_00 + CF_03) + ((3.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_0 +
            (((1.0/4.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_1 +
            Recon_0;

       beta_0 = ((1.0/4.0))*((CS_12*max_lambda_11 + CF_12)*(CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (CS_12*max_lambda_11 + CF_12))*(CS_11*max_lambda_11
            + CF_11);

       beta_1 = ((1.0/4.0))*((CS_11*max_lambda_11 + CF_11)*(CS_11*max_lambda_11 + CF_11)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - (CS_11*max_lambda_11 + CF_11))*(CS_10*max_lambda_11
            + CF_10);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_1 = (-(1.0/4.0)*(CS_10*max_lambda_11 + CF_10) + ((3.0/4.0))*(CS_11*max_lambda_11 + CF_11))*omega_1 +
            (((1.0/4.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(CS_12*max_lambda_11 + CF_12))*omega_0 + Recon_1;

       beta_0 = ((1.0/4.0))*((-CS_13*max_lambda_11 + CF_13)*(-CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (-CS_13*max_lambda_11 +
            CF_13))*(-CS_12*max_lambda_11 + CF_12);

       beta_1 = ((1.0/4.0))*((-CS_12*max_lambda_11 + CF_12)*(-CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - (-CS_12*max_lambda_11 +
            CF_12))*(-CS_11*max_lambda_11 + CF_11);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = fmax(rj1, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_1 = (-(1.0/4.0)*(-CS_13*max_lambda_11 + CF_13) + ((3.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_0 +
            (((1.0/4.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_1 +
            Recon_1;

       beta_0 = ((1.0/4.0))*((CS_22*max_lambda_22 + CF_22)*(CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (CS_22*max_lambda_22 + CF_22))*(CS_21*max_lambda_22
            + CF_21);

       beta_1 = ((1.0/4.0))*((CS_21*max_lambda_22 + CF_21)*(CS_21*max_lambda_22 + CF_21)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - (CS_21*max_lambda_22 + CF_21))*(CS_20*max_lambda_22
            + CF_20);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_2 = (-(1.0/4.0)*(CS_20*max_lambda_22 + CF_20) + ((3.0/4.0))*(CS_21*max_lambda_22 + CF_21))*omega_1 +
            (((1.0/4.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(CS_22*max_lambda_22 + CF_22))*omega_0 + Recon_2;

       beta_0 = ((1.0/4.0))*((-CS_23*max_lambda_22 + CF_23)*(-CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (-CS_23*max_lambda_22 +
            CF_23))*(-CS_22*max_lambda_22 + CF_22);

       beta_1 = ((1.0/4.0))*((-CS_22*max_lambda_22 + CF_22)*(-CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - (-CS_22*max_lambda_22 +
            CF_22))*(-CS_21*max_lambda_22 + CF_21);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = fmax(rj2, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_2 = (-(1.0/4.0)*(-CS_23*max_lambda_22 + CF_23) + ((3.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_0 +
            (((1.0/4.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_1 +
            Recon_2;

       beta_0 = ((1.0/4.0))*((CS_32*max_lambda_33 + CF_32)*(CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (CS_32*max_lambda_33 + CF_32))*(CS_31*max_lambda_33
            + CF_31);

       beta_1 = ((1.0/4.0))*((CS_31*max_lambda_33 + CF_31)*(CS_31*max_lambda_33 + CF_31)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - (CS_31*max_lambda_33 + CF_31))*(CS_30*max_lambda_33
            + CF_30);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_3 = (-(1.0/4.0)*(CS_30*max_lambda_33 + CF_30) + ((3.0/4.0))*(CS_31*max_lambda_33 + CF_31))*omega_1 +
            (((1.0/4.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(CS_32*max_lambda_33 + CF_32))*omega_0 + Recon_3;

       beta_0 = ((1.0/4.0))*((-CS_33*max_lambda_33 + CF_33)*(-CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (-CS_33*max_lambda_33 +
            CF_33))*(-CS_32*max_lambda_33 + CF_32);

       beta_1 = ((1.0/4.0))*((-CS_32*max_lambda_33 + CF_32)*(-CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - (-CS_32*max_lambda_33 +
            CF_32))*(-CS_31*max_lambda_33 + CF_31);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = fmax(rj3, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_3 = (-(1.0/4.0)*(-CS_33*max_lambda_33 + CF_33) + ((3.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_0 +
            (((1.0/4.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_1 +
            Recon_3;

       beta_0 = ((1.0/4.0))*((CS_42*max_lambda_44 + CF_42)*(CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (CS_42*max_lambda_44 + CF_42))*(CS_41*max_lambda_44
            + CF_41);

       beta_1 = ((1.0/4.0))*((CS_41*max_lambda_44 + CF_41)*(CS_41*max_lambda_44 + CF_41)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - (CS_41*max_lambda_44 + CF_41))*(CS_40*max_lambda_44
            + CF_40);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_4 = (-(1.0/4.0)*(CS_40*max_lambda_44 + CF_40) + ((3.0/4.0))*(CS_41*max_lambda_44 + CF_41))*omega_1 +
            (((1.0/4.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(CS_42*max_lambda_44 + CF_42))*omega_0 + Recon_4;

       beta_0 = ((1.0/4.0))*((-CS_43*max_lambda_44 + CF_43)*(-CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (-CS_43*max_lambda_44 +
            CF_43))*(-CS_42*max_lambda_44 + CF_42);

       beta_1 = ((1.0/4.0))*((-CS_42*max_lambda_44 + CF_42)*(-CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - (-CS_42*max_lambda_44 +
            CF_42))*(-CS_41*max_lambda_44 + CF_41);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = fmax(rj4, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_4 = (-(1.0/4.0)*(-CS_43*max_lambda_44 + CF_43) + ((3.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_0 +
            (((1.0/4.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_1 +
            Recon_4;

      Recon_0 = (((1.0/12.0))*(-7*CF_01 - 7*CF_02 + CF_00 + CF_03) + Recon_0)*rj0;

      Recon_1 = (((1.0/12.0))*(-7*CF_11 - 7*CF_12 + CF_10 + CF_13) + Recon_1)*rj1;

      Recon_2 = (((1.0/12.0))*(-7*CF_21 - 7*CF_22 + CF_20 + CF_23) + Recon_2)*rj2;

      Recon_3 = (((1.0/12.0))*(-7*CF_31 - 7*CF_32 + CF_30 + CF_33) + Recon_3)*rj3;

      Recon_4 = (((1.0/12.0))*(-7*CF_41 - 7*CF_42 + CF_40 + CF_43) + Recon_4)*rj4;

       wk0_B0(0,0,0) = 0.707106781186547*AVG_0_rho*Recon_3*inv_AVG_a + 0.707106781186547*AVG_0_rho*Recon_4*inv_AVG_a +
            Recon_0;

       wk1_B0(0,0,0) = AVG_0_u0*Recon_0 + 0.707106781186547*(-AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_0_a + AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a;

       wk2_B0(0,0,0) = AVG_0_u1*Recon_0 - AVG_0_rho*Recon_2 + 0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_0_rho*AVG_0_u1*Recon_4*inv_AVG_a;

       wk3_B0(0,0,0) = AVG_0_rho*Recon_1 + AVG_0_u2*Recon_0 + 0.707106781186547*AVG_0_rho*AVG_0_u2*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_0_rho*AVG_0_u2*Recon_4*inv_AVG_a;

       wk4_B0(0,0,0) = (((1.0/2.0))*(AVG_0_u0*AVG_0_u0) + ((1.0/2.0))*(AVG_0_u1*AVG_0_u1) +
            ((1.0/2.0))*(AVG_0_u2*AVG_0_u2))*Recon_0 + AVG_0_rho*AVG_0_u2*Recon_1 - AVG_0_rho*AVG_0_u1*Recon_2 +
            0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) +
            (AVG_0_u2*AVG_0_u2))*gamma_m1)*invgamma_m1 + AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_0_a*AVG_0_a) + ((1.0/2.0))*((AVG_0_u0*AVG_0_u0) + (AVG_0_u1*AVG_0_u1) +
            (AVG_0_u2*AVG_0_u2))*gamma_m1)*invgamma_m1 - AVG_0_a*AVG_0_u0)*AVG_0_rho*Recon_4*inv_AVG_a;

   }

   else{

      wk0_B0(0,0,0) = 0.0;

      wk1_B0(0,0,0) = 0.0;

      wk2_B0(0,0,0) = 0.0;

      wk3_B0(0,0,0) = 0.0;

      wk4_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel047(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&Residual0_B0, ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double> &Residual3_B0, ACC<double>
&Residual4_B0)
{
   double AVG_1_1_LEV_00 = 0.0;
   double AVG_1_1_LEV_03 = 0.0;
   double AVG_1_1_LEV_10 = 0.0;
   double AVG_1_1_LEV_11 = 0.0;
   double AVG_1_1_LEV_12 = 0.0;
   double AVG_1_1_LEV_13 = 0.0;
   double AVG_1_1_LEV_14 = 0.0;
   double AVG_1_1_LEV_20 = 0.0;
   double AVG_1_1_LEV_21 = 0.0;
   double AVG_1_1_LEV_30 = 0.0;
   double AVG_1_1_LEV_31 = 0.0;
   double AVG_1_1_LEV_32 = 0.0;
   double AVG_1_1_LEV_33 = 0.0;
   double AVG_1_1_LEV_34 = 0.0;
   double AVG_1_1_LEV_40 = 0.0;
   double AVG_1_1_LEV_41 = 0.0;
   double AVG_1_1_LEV_42 = 0.0;
   double AVG_1_1_LEV_43 = 0.0;
   double AVG_1_1_LEV_44 = 0.0;
   double AVG_1_a = 0.0;
   double AVG_1_inv_rho = 0.0;
   double AVG_1_rho = 0.0;
   double AVG_1_u0 = 0.0;
   double AVG_1_u1 = 0.0;
   double AVG_1_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(0,-2,0), fmax(kappa_B0(0,-1,0), fmax(kappa_B0(0,0,0), fmax(kappa_B0(0,1,0), fmax(kappa_B0(0,2,0),
      kappa_B0(0,-3,0)))))) > Ducros_check){

      AVG_1_rho = sqrt((rho_B0(0,0,0)*rho_B0(0,1,0)));

      AVG_1_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(0,1,0))));

      AVG_1_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u0_B0(0,1,0))*AVG_1_inv_rho;

      AVG_1_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u1_B0(0,1,0))*AVG_1_inv_rho;

      AVG_1_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(0,1,0))*u2_B0(0,1,0))*AVG_1_inv_rho;

       AVG_1_a = sqrt(((-(1.0/2.0)*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) + (AVG_1_u2*AVG_1_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(0,1,0) +
            rhoE_B0(0,1,0))/sqrt(rho_B0(0,1,0)))*AVG_1_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_1_a);

      inv_AVG_rho = 1.0/(AVG_1_rho);

      AVG_1_1_LEV_00 = AVG_1_u2*inv_AVG_rho;

      AVG_1_1_LEV_03 = -inv_AVG_rho;

       AVG_1_1_LEV_10 = -(1.0/2.0)*(-2 - (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_1_u2*AVG_1_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_1_u0*AVG_1_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_1_u1*AVG_1_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_1_u2*AVG_1_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_1_1_LEV_11 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u0;

      AVG_1_1_LEV_12 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u1;

      AVG_1_1_LEV_13 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_1_u2;

      AVG_1_1_LEV_14 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

      AVG_1_1_LEV_20 = -AVG_1_u0*inv_AVG_rho;

      AVG_1_1_LEV_21 = inv_AVG_rho;

       AVG_1_1_LEV_30 = -0.353553390593274*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) + (AVG_1_u2*AVG_1_u2) -
            (AVG_1_u0*AVG_1_u0)*gama - (AVG_1_u1*AVG_1_u1)*gama - (AVG_1_u2*AVG_1_u2)*gama +
            2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_31 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_32 = 0.707106781186547*(-gama*AVG_1_u1 + AVG_1_a + AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_33 = -0.707106781186547*gamma_m1*AVG_1_u2*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_1_1_LEV_40 = 0.353553390593274*(-(AVG_1_u0*AVG_1_u0) - (AVG_1_u1*AVG_1_u1) - (AVG_1_u2*AVG_1_u2) +
            (AVG_1_u0*AVG_1_u0)*gama + (AVG_1_u1*AVG_1_u1)*gama + (AVG_1_u2*AVG_1_u2)*gama +
            2*AVG_1_a*AVG_1_u1)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_41 = -0.707106781186547*gamma_m1*AVG_1_u0*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_42 = -0.707106781186547*(-AVG_1_u1 + gama*AVG_1_u1 + AVG_1_a)*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_43 = -0.707106781186547*gamma_m1*AVG_1_u2*inv_AVG_a*inv_AVG_rho;

      AVG_1_1_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

      CF_00 = rhou1_B0(0,-1,0)*AVG_1_1_LEV_00 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_03;

       CF_10 = p_B0(0,-1,0)*AVG_1_1_LEV_12 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_10 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_14
            + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 +
            u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

      CF_20 = rhou1_B0(0,-1,0)*AVG_1_1_LEV_20 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_21;

       CF_30 = p_B0(0,-1,0)*AVG_1_1_LEV_32 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_30 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_34
            + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 +
            u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

       CF_40 = p_B0(0,-1,0)*AVG_1_1_LEV_42 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_40 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_44
            + u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 +
            u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

      CS_00 = rho_B0(0,-1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_03;

       CS_10 = rho_B0(0,-1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

      CS_20 = rho_B0(0,-1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_21;

       CS_30 = rho_B0(0,-1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

       CS_40 = rho_B0(0,-1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

      CF_01 = rhou1_B0(0,0,0)*AVG_1_1_LEV_00 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_03;

       CF_11 = p_B0(0,0,0)*AVG_1_1_LEV_12 + rhou1_B0(0,0,0)*AVG_1_1_LEV_10 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_14 +
            u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_14 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_11 +
            u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_12 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

      CF_21 = rhou1_B0(0,0,0)*AVG_1_1_LEV_20 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_21;

       CF_31 = p_B0(0,0,0)*AVG_1_1_LEV_32 + rhou1_B0(0,0,0)*AVG_1_1_LEV_30 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_34 +
            u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_34 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_31 +
            u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_32 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

       CF_41 = p_B0(0,0,0)*AVG_1_1_LEV_42 + rhou1_B0(0,0,0)*AVG_1_1_LEV_40 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_44 +
            u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_44 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_41 +
            u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_42 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

      CS_01 = rho_B0(0,0,0)*AVG_1_1_LEV_00 + rhou2_B0(0,0,0)*AVG_1_1_LEV_03;

       CS_11 = rho_B0(0,0,0)*AVG_1_1_LEV_10 + rhoE_B0(0,0,0)*AVG_1_1_LEV_14 + rhou0_B0(0,0,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_12 + rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

      CS_21 = rho_B0(0,0,0)*AVG_1_1_LEV_20 + rhou0_B0(0,0,0)*AVG_1_1_LEV_21;

       CS_31 = rho_B0(0,0,0)*AVG_1_1_LEV_30 + rhoE_B0(0,0,0)*AVG_1_1_LEV_34 + rhou0_B0(0,0,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_32 + rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

       CS_41 = rho_B0(0,0,0)*AVG_1_1_LEV_40 + rhoE_B0(0,0,0)*AVG_1_1_LEV_44 + rhou0_B0(0,0,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,0,0)*AVG_1_1_LEV_42 + rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

      CF_02 = rhou1_B0(0,1,0)*AVG_1_1_LEV_00 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_03;

       CF_12 = p_B0(0,1,0)*AVG_1_1_LEV_12 + rhou1_B0(0,1,0)*AVG_1_1_LEV_10 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_14 +
            u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_14 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_11 +
            u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_12 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

      CF_22 = rhou1_B0(0,1,0)*AVG_1_1_LEV_20 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_21;

       CF_32 = p_B0(0,1,0)*AVG_1_1_LEV_32 + rhou1_B0(0,1,0)*AVG_1_1_LEV_30 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_34 +
            u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_34 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_31 +
            u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_32 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

       CF_42 = p_B0(0,1,0)*AVG_1_1_LEV_42 + rhou1_B0(0,1,0)*AVG_1_1_LEV_40 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_44 +
            u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_44 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_41 +
            u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_42 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

      CS_02 = rho_B0(0,1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,1,0)*AVG_1_1_LEV_03;

       CS_12 = rho_B0(0,1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,1,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

      CS_22 = rho_B0(0,1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,1,0)*AVG_1_1_LEV_21;

       CS_32 = rho_B0(0,1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,1,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

       CS_42 = rho_B0(0,1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,1,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

      CF_03 = rhou1_B0(0,2,0)*AVG_1_1_LEV_00 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_03;

       CF_13 = p_B0(0,2,0)*AVG_1_1_LEV_12 + rhou1_B0(0,2,0)*AVG_1_1_LEV_10 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_14 +
            u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_14 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_11 +
            u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_12 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

      CF_23 = rhou1_B0(0,2,0)*AVG_1_1_LEV_20 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_21;

       CF_33 = p_B0(0,2,0)*AVG_1_1_LEV_32 + rhou1_B0(0,2,0)*AVG_1_1_LEV_30 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_34 +
            u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_34 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_31 +
            u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_32 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

       CF_43 = p_B0(0,2,0)*AVG_1_1_LEV_42 + rhou1_B0(0,2,0)*AVG_1_1_LEV_40 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_44 +
            u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_44 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_41 +
            u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_42 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

      CS_03 = rho_B0(0,2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,2,0)*AVG_1_1_LEV_03;

       CS_13 = rho_B0(0,2,0)*AVG_1_1_LEV_10 + rhoE_B0(0,2,0)*AVG_1_1_LEV_14 + rhou0_B0(0,2,0)*AVG_1_1_LEV_11 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_12 + rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

      CS_23 = rho_B0(0,2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,2,0)*AVG_1_1_LEV_21;

       CS_33 = rho_B0(0,2,0)*AVG_1_1_LEV_30 + rhoE_B0(0,2,0)*AVG_1_1_LEV_34 + rhou0_B0(0,2,0)*AVG_1_1_LEV_31 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_32 + rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

       CS_43 = rho_B0(0,2,0)*AVG_1_1_LEV_40 + rhoE_B0(0,2,0)*AVG_1_1_LEV_44 + rhou0_B0(0,2,0)*AVG_1_1_LEV_41 +
            rhou1_B0(0,2,0)*AVG_1_1_LEV_42 + rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u1_B0(0,0,0)), fabs(u1_B0(0,1,0)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,0,0) + u1_B0(0,0,0)), fabs(a_B0(0,1,0) + u1_B0(0,1,0)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u1_B0(0,0,0) + a_B0(0,0,0)), fabs(-u1_B0(0,1,0) + a_B0(0,1,0)));

       beta_0 = ((1.0/4.0))*((CS_02*max_lambda_00 + CF_02)*(CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (CS_02*max_lambda_00 + CF_02))*(CS_01*max_lambda_00
            + CF_01);

       beta_1 = ((1.0/4.0))*((CS_01*max_lambda_00 + CF_01)*(CS_01*max_lambda_00 + CF_01)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - (CS_01*max_lambda_00 + CF_01))*(CS_00*max_lambda_00
            + CF_00);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_0 = (-(1.0/4.0)*(CS_00*max_lambda_00 + CF_00) + ((3.0/4.0))*(CS_01*max_lambda_00 + CF_01))*omega_1 +
            (((1.0/4.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(CS_02*max_lambda_00 + CF_02))*omega_0 + Recon_0;

       beta_0 = ((1.0/4.0))*((-CS_03*max_lambda_00 + CF_03)*(-CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (-CS_03*max_lambda_00 +
            CF_03))*(-CS_02*max_lambda_00 + CF_02);

       beta_1 = ((1.0/4.0))*((-CS_02*max_lambda_00 + CF_02)*(-CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - (-CS_02*max_lambda_00 +
            CF_02))*(-CS_01*max_lambda_00 + CF_01);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = fmax(rj0, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_0 = (-(1.0/4.0)*(-CS_03*max_lambda_00 + CF_03) + ((3.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_0 +
            (((1.0/4.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_1 +
            Recon_0;

       beta_0 = ((1.0/4.0))*((CS_12*max_lambda_11 + CF_12)*(CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (CS_12*max_lambda_11 + CF_12))*(CS_11*max_lambda_11
            + CF_11);

       beta_1 = ((1.0/4.0))*((CS_11*max_lambda_11 + CF_11)*(CS_11*max_lambda_11 + CF_11)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - (CS_11*max_lambda_11 + CF_11))*(CS_10*max_lambda_11
            + CF_10);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_1 = (-(1.0/4.0)*(CS_10*max_lambda_11 + CF_10) + ((3.0/4.0))*(CS_11*max_lambda_11 + CF_11))*omega_1 +
            (((1.0/4.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(CS_12*max_lambda_11 + CF_12))*omega_0 + Recon_1;

       beta_0 = ((1.0/4.0))*((-CS_13*max_lambda_11 + CF_13)*(-CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (-CS_13*max_lambda_11 +
            CF_13))*(-CS_12*max_lambda_11 + CF_12);

       beta_1 = ((1.0/4.0))*((-CS_12*max_lambda_11 + CF_12)*(-CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - (-CS_12*max_lambda_11 +
            CF_12))*(-CS_11*max_lambda_11 + CF_11);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = fmax(rj1, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_1 = (-(1.0/4.0)*(-CS_13*max_lambda_11 + CF_13) + ((3.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_0 +
            (((1.0/4.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_1 +
            Recon_1;

       beta_0 = ((1.0/4.0))*((CS_22*max_lambda_22 + CF_22)*(CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (CS_22*max_lambda_22 + CF_22))*(CS_21*max_lambda_22
            + CF_21);

       beta_1 = ((1.0/4.0))*((CS_21*max_lambda_22 + CF_21)*(CS_21*max_lambda_22 + CF_21)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - (CS_21*max_lambda_22 + CF_21))*(CS_20*max_lambda_22
            + CF_20);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_2 = (-(1.0/4.0)*(CS_20*max_lambda_22 + CF_20) + ((3.0/4.0))*(CS_21*max_lambda_22 + CF_21))*omega_1 +
            (((1.0/4.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(CS_22*max_lambda_22 + CF_22))*omega_0 + Recon_2;

       beta_0 = ((1.0/4.0))*((-CS_23*max_lambda_22 + CF_23)*(-CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (-CS_23*max_lambda_22 +
            CF_23))*(-CS_22*max_lambda_22 + CF_22);

       beta_1 = ((1.0/4.0))*((-CS_22*max_lambda_22 + CF_22)*(-CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - (-CS_22*max_lambda_22 +
            CF_22))*(-CS_21*max_lambda_22 + CF_21);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = fmax(rj2, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_2 = (-(1.0/4.0)*(-CS_23*max_lambda_22 + CF_23) + ((3.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_0 +
            (((1.0/4.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_1 +
            Recon_2;

       beta_0 = ((1.0/4.0))*((CS_32*max_lambda_33 + CF_32)*(CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (CS_32*max_lambda_33 + CF_32))*(CS_31*max_lambda_33
            + CF_31);

       beta_1 = ((1.0/4.0))*((CS_31*max_lambda_33 + CF_31)*(CS_31*max_lambda_33 + CF_31)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - (CS_31*max_lambda_33 + CF_31))*(CS_30*max_lambda_33
            + CF_30);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_3 = (-(1.0/4.0)*(CS_30*max_lambda_33 + CF_30) + ((3.0/4.0))*(CS_31*max_lambda_33 + CF_31))*omega_1 +
            (((1.0/4.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(CS_32*max_lambda_33 + CF_32))*omega_0 + Recon_3;

       beta_0 = ((1.0/4.0))*((-CS_33*max_lambda_33 + CF_33)*(-CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (-CS_33*max_lambda_33 +
            CF_33))*(-CS_32*max_lambda_33 + CF_32);

       beta_1 = ((1.0/4.0))*((-CS_32*max_lambda_33 + CF_32)*(-CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - (-CS_32*max_lambda_33 +
            CF_32))*(-CS_31*max_lambda_33 + CF_31);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = fmax(rj3, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_3 = (-(1.0/4.0)*(-CS_33*max_lambda_33 + CF_33) + ((3.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_0 +
            (((1.0/4.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_1 +
            Recon_3;

       beta_0 = ((1.0/4.0))*((CS_42*max_lambda_44 + CF_42)*(CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (CS_42*max_lambda_44 + CF_42))*(CS_41*max_lambda_44
            + CF_41);

       beta_1 = ((1.0/4.0))*((CS_41*max_lambda_44 + CF_41)*(CS_41*max_lambda_44 + CF_41)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - (CS_41*max_lambda_44 + CF_41))*(CS_40*max_lambda_44
            + CF_40);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_4 = (-(1.0/4.0)*(CS_40*max_lambda_44 + CF_40) + ((3.0/4.0))*(CS_41*max_lambda_44 + CF_41))*omega_1 +
            (((1.0/4.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(CS_42*max_lambda_44 + CF_42))*omega_0 + Recon_4;

       beta_0 = ((1.0/4.0))*((-CS_43*max_lambda_44 + CF_43)*(-CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (-CS_43*max_lambda_44 +
            CF_43))*(-CS_42*max_lambda_44 + CF_42);

       beta_1 = ((1.0/4.0))*((-CS_42*max_lambda_44 + CF_42)*(-CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - (-CS_42*max_lambda_44 +
            CF_42))*(-CS_41*max_lambda_44 + CF_41);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = fmax(rj4, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_4 = (-(1.0/4.0)*(-CS_43*max_lambda_44 + CF_43) + ((3.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_0 +
            (((1.0/4.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_1 +
            Recon_4;

      Recon_0 = (((1.0/12.0))*(-7*CF_01 - 7*CF_02 + CF_00 + CF_03) + Recon_0)*rj0;

      Recon_1 = (((1.0/12.0))*(-7*CF_11 - 7*CF_12 + CF_10 + CF_13) + Recon_1)*rj1;

      Recon_2 = (((1.0/12.0))*(-7*CF_21 - 7*CF_22 + CF_20 + CF_23) + Recon_2)*rj2;

      Recon_3 = (((1.0/12.0))*(-7*CF_31 - 7*CF_32 + CF_30 + CF_33) + Recon_3)*rj3;

      Recon_4 = (((1.0/12.0))*(-7*CF_41 - 7*CF_42 + CF_40 + CF_43) + Recon_4)*rj4;

       Residual0_B0(0,0,0) = 0.707106781186547*AVG_1_rho*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*Recon_4*inv_AVG_a + Recon_1;

       Residual1_B0(0,0,0) = AVG_1_rho*Recon_2 + AVG_1_u0*Recon_1 +
            0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_4*inv_AVG_a;

       Residual2_B0(0,0,0) = AVG_1_u1*Recon_1 + 0.707106781186547*(-AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a;

       Residual3_B0(0,0,0) = AVG_1_u2*Recon_1 - AVG_1_rho*Recon_0 +
            0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_4*inv_AVG_a;

       Residual4_B0(0,0,0) = (((1.0/2.0))*(AVG_1_u0*AVG_1_u0) + ((1.0/2.0))*(AVG_1_u1*AVG_1_u1) +
            ((1.0/2.0))*(AVG_1_u2*AVG_1_u2))*Recon_1 + AVG_1_rho*AVG_1_u0*Recon_2 - AVG_1_rho*AVG_1_u2*Recon_0 +
            0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
            (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 + AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
            (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 - AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a;

   }

   else{

      Residual0_B0(0,0,0) = 0.0;

      Residual1_B0(0,0,0) = 0.0;

      Residual2_B0(0,0,0) = 0.0;

      Residual3_B0(0,0,0) = 0.0;

      Residual4_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel048(const ACC<double> &a_B0, const ACC<double> &kappa_B0, const ACC<double> &p_B0, const
ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const
ACC<double> &rhou2_B0, const ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double>
&rhoE_RKold_B0, ACC<double> &rho_RKold_B0, ACC<double> &rhou0_RKold_B0, ACC<double> &rhou1_RKold_B0, ACC<double>
&rhou2_RKold_B0)
{
   double AVG_2_2_LEV_00 = 0.0;
   double AVG_2_2_LEV_02 = 0.0;
   double AVG_2_2_LEV_10 = 0.0;
   double AVG_2_2_LEV_11 = 0.0;
   double AVG_2_2_LEV_20 = 0.0;
   double AVG_2_2_LEV_21 = 0.0;
   double AVG_2_2_LEV_22 = 0.0;
   double AVG_2_2_LEV_23 = 0.0;
   double AVG_2_2_LEV_24 = 0.0;
   double AVG_2_2_LEV_30 = 0.0;
   double AVG_2_2_LEV_31 = 0.0;
   double AVG_2_2_LEV_32 = 0.0;
   double AVG_2_2_LEV_33 = 0.0;
   double AVG_2_2_LEV_34 = 0.0;
   double AVG_2_2_LEV_40 = 0.0;
   double AVG_2_2_LEV_41 = 0.0;
   double AVG_2_2_LEV_42 = 0.0;
   double AVG_2_2_LEV_43 = 0.0;
   double AVG_2_2_LEV_44 = 0.0;
   double AVG_2_a = 0.0;
   double AVG_2_inv_rho = 0.0;
   double AVG_2_rho = 0.0;
   double AVG_2_u0 = 0.0;
   double AVG_2_u1 = 0.0;
   double AVG_2_u2 = 0.0;
   double CF_00 = 0.0;
   double CF_01 = 0.0;
   double CF_02 = 0.0;
   double CF_03 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
   double omega_0 = 0.0;
   double omega_1 = 0.0;
   double rj0 = 0.0;
   double rj1 = 0.0;
   double rj2 = 0.0;
   double rj3 = 0.0;
   double rj4 = 0.0;
    if (fmax(kappa_B0(0,0,1), fmax(kappa_B0(0,0,0), fmax(kappa_B0(0,0,2), fmax(kappa_B0(0,0,-3), fmax(kappa_B0(0,0,-2),
      kappa_B0(0,0,-1)))))) > Ducros_check){

      AVG_2_rho = sqrt((rho_B0(0,0,0)*rho_B0(0,0,1)));

      AVG_2_inv_rho = 1.0/((sqrt(rho_B0(0,0,0)) + sqrt(rho_B0(0,0,1))));

      AVG_2_u0 = (sqrt(rho_B0(0,0,0))*u0_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u0_B0(0,0,1))*AVG_2_inv_rho;

      AVG_2_u1 = (sqrt(rho_B0(0,0,0))*u1_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u1_B0(0,0,1))*AVG_2_inv_rho;

      AVG_2_u2 = (sqrt(rho_B0(0,0,0))*u2_B0(0,0,0) + sqrt(rho_B0(0,0,1))*u2_B0(0,0,1))*AVG_2_inv_rho;

       AVG_2_a = sqrt(((-(1.0/2.0)*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) + (AVG_2_u2*AVG_2_u2)) + ((p_B0(0,0,0) +
            rhoE_B0(0,0,0))/sqrt(rho_B0(0,0,0)) + (p_B0(0,0,1) +
            rhoE_B0(0,0,1))/sqrt(rho_B0(0,0,1)))*AVG_2_inv_rho)*gamma_m1));

      inv_AVG_a = 1.0/(AVG_2_a);

      inv_AVG_rho = 1.0/(AVG_2_rho);

      AVG_2_2_LEV_00 = -AVG_2_u1*inv_AVG_rho;

      AVG_2_2_LEV_02 = inv_AVG_rho;

      AVG_2_2_LEV_10 = AVG_2_u0*inv_AVG_rho;

      AVG_2_2_LEV_11 = -inv_AVG_rho;

       AVG_2_2_LEV_20 = -(1.0/2.0)*(-2 - (AVG_2_u0*AVG_2_u0)*(inv_AVG_a*inv_AVG_a) -
            (AVG_2_u1*AVG_2_u1)*(inv_AVG_a*inv_AVG_a) - (AVG_2_u2*AVG_2_u2)*(inv_AVG_a*inv_AVG_a) +
            (AVG_2_u0*AVG_2_u0)*(inv_AVG_a*inv_AVG_a)*gama + (AVG_2_u1*AVG_2_u1)*(inv_AVG_a*inv_AVG_a)*gama +
            (AVG_2_u2*AVG_2_u2)*(inv_AVG_a*inv_AVG_a)*gama);

      AVG_2_2_LEV_21 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u0;

      AVG_2_2_LEV_22 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u1;

      AVG_2_2_LEV_23 = (inv_AVG_a*inv_AVG_a)*gamma_m1*AVG_2_u2;

      AVG_2_2_LEV_24 = -(inv_AVG_a*inv_AVG_a)*gamma_m1;

       AVG_2_2_LEV_30 = -0.353553390593274*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) + (AVG_2_u2*AVG_2_u2) -
            (AVG_2_u0*AVG_2_u0)*gama - (AVG_2_u1*AVG_2_u1)*gama - (AVG_2_u2*AVG_2_u2)*gama +
            2*AVG_2_a*AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_31 = -0.707106781186547*gamma_m1*AVG_2_u0*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_32 = -0.707106781186547*gamma_m1*AVG_2_u1*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_33 = 0.707106781186547*(-gama*AVG_2_u2 + AVG_2_a + AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_34 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

       AVG_2_2_LEV_40 = 0.353553390593274*(-(AVG_2_u0*AVG_2_u0) - (AVG_2_u1*AVG_2_u1) - (AVG_2_u2*AVG_2_u2) +
            (AVG_2_u0*AVG_2_u0)*gama + (AVG_2_u1*AVG_2_u1)*gama + (AVG_2_u2*AVG_2_u2)*gama +
            2*AVG_2_a*AVG_2_u2)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_41 = -0.707106781186547*gamma_m1*AVG_2_u0*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_42 = -0.707106781186547*gamma_m1*AVG_2_u1*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_43 = -0.707106781186547*(-AVG_2_u2 + gama*AVG_2_u2 + AVG_2_a)*inv_AVG_a*inv_AVG_rho;

      AVG_2_2_LEV_44 = 0.707106781186547*gamma_m1*inv_AVG_a*inv_AVG_rho;

      CF_00 = rhou2_B0(0,0,-1)*AVG_2_2_LEV_00 + u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_02;

      CF_10 = rhou2_B0(0,0,-1)*AVG_2_2_LEV_10 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_11;

       CF_20 = p_B0(0,0,-1)*AVG_2_2_LEV_23 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_20 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_24
            + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 +
            u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

       CF_30 = p_B0(0,0,-1)*AVG_2_2_LEV_33 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_30 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_34
            + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 +
            u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

       CF_40 = p_B0(0,0,-1)*AVG_2_2_LEV_43 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_40 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_44
            + u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 +
            u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

      CS_00 = rho_B0(0,0,-1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-1)*AVG_2_2_LEV_02;

      CS_10 = rho_B0(0,0,-1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_11;

       CS_20 = rho_B0(0,0,-1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

       CS_30 = rho_B0(0,0,-1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

       CS_40 = rho_B0(0,0,-1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

      CF_01 = rhou2_B0(0,0,0)*AVG_2_2_LEV_00 + u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_02;

      CF_11 = rhou2_B0(0,0,0)*AVG_2_2_LEV_10 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_11;

       CF_21 = p_B0(0,0,0)*AVG_2_2_LEV_23 + rhou2_B0(0,0,0)*AVG_2_2_LEV_20 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_24 +
            u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_24 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_21 +
            u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_22 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

       CF_31 = p_B0(0,0,0)*AVG_2_2_LEV_33 + rhou2_B0(0,0,0)*AVG_2_2_LEV_30 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_34 +
            u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_34 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_31 +
            u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_32 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

       CF_41 = p_B0(0,0,0)*AVG_2_2_LEV_43 + rhou2_B0(0,0,0)*AVG_2_2_LEV_40 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_44 +
            u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_44 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_41 +
            u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_42 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

      CS_01 = rho_B0(0,0,0)*AVG_2_2_LEV_00 + rhou1_B0(0,0,0)*AVG_2_2_LEV_02;

      CS_11 = rho_B0(0,0,0)*AVG_2_2_LEV_10 + rhou0_B0(0,0,0)*AVG_2_2_LEV_11;

       CS_21 = rho_B0(0,0,0)*AVG_2_2_LEV_20 + rhoE_B0(0,0,0)*AVG_2_2_LEV_24 + rhou0_B0(0,0,0)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_22 + rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

       CS_31 = rho_B0(0,0,0)*AVG_2_2_LEV_30 + rhoE_B0(0,0,0)*AVG_2_2_LEV_34 + rhou0_B0(0,0,0)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_32 + rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

       CS_41 = rho_B0(0,0,0)*AVG_2_2_LEV_40 + rhoE_B0(0,0,0)*AVG_2_2_LEV_44 + rhou0_B0(0,0,0)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,0)*AVG_2_2_LEV_42 + rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

      CF_02 = rhou2_B0(0,0,1)*AVG_2_2_LEV_00 + u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_02;

      CF_12 = rhou2_B0(0,0,1)*AVG_2_2_LEV_10 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_11;

       CF_22 = p_B0(0,0,1)*AVG_2_2_LEV_23 + rhou2_B0(0,0,1)*AVG_2_2_LEV_20 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_24 +
            u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_24 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_21 +
            u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_22 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

       CF_32 = p_B0(0,0,1)*AVG_2_2_LEV_33 + rhou2_B0(0,0,1)*AVG_2_2_LEV_30 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_34 +
            u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_34 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_31 +
            u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_32 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

       CF_42 = p_B0(0,0,1)*AVG_2_2_LEV_43 + rhou2_B0(0,0,1)*AVG_2_2_LEV_40 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_44 +
            u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_44 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_41 +
            u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_42 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

      CS_02 = rho_B0(0,0,1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,1)*AVG_2_2_LEV_02;

      CS_12 = rho_B0(0,0,1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,1)*AVG_2_2_LEV_11;

       CS_22 = rho_B0(0,0,1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,1)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

       CS_32 = rho_B0(0,0,1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,1)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

       CS_42 = rho_B0(0,0,1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,1)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

      CF_03 = rhou2_B0(0,0,2)*AVG_2_2_LEV_00 + u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_02;

      CF_13 = rhou2_B0(0,0,2)*AVG_2_2_LEV_10 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_11;

       CF_23 = p_B0(0,0,2)*AVG_2_2_LEV_23 + rhou2_B0(0,0,2)*AVG_2_2_LEV_20 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_24 +
            u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_24 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_21 +
            u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_22 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

       CF_33 = p_B0(0,0,2)*AVG_2_2_LEV_33 + rhou2_B0(0,0,2)*AVG_2_2_LEV_30 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_34 +
            u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_34 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_31 +
            u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_32 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

       CF_43 = p_B0(0,0,2)*AVG_2_2_LEV_43 + rhou2_B0(0,0,2)*AVG_2_2_LEV_40 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_44 +
            u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_44 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_41 +
            u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_42 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

      CS_03 = rho_B0(0,0,2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,2)*AVG_2_2_LEV_02;

      CS_13 = rho_B0(0,0,2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,2)*AVG_2_2_LEV_11;

       CS_23 = rho_B0(0,0,2)*AVG_2_2_LEV_20 + rhoE_B0(0,0,2)*AVG_2_2_LEV_24 + rhou0_B0(0,0,2)*AVG_2_2_LEV_21 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_22 + rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

       CS_33 = rho_B0(0,0,2)*AVG_2_2_LEV_30 + rhoE_B0(0,0,2)*AVG_2_2_LEV_34 + rhou0_B0(0,0,2)*AVG_2_2_LEV_31 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_32 + rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

       CS_43 = rho_B0(0,0,2)*AVG_2_2_LEV_40 + rhoE_B0(0,0,2)*AVG_2_2_LEV_44 + rhou0_B0(0,0,2)*AVG_2_2_LEV_41 +
            rhou1_B0(0,0,2)*AVG_2_2_LEV_42 + rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

      max_lambda_00 = shock_filter_control*fmax(fabs(u2_B0(0,0,0)), fabs(u2_B0(0,0,1)));

      max_lambda_11 = max_lambda_00;

      max_lambda_22 = max_lambda_00;

      max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,0,0) + u2_B0(0,0,0)), fabs(a_B0(0,0,1) + u2_B0(0,0,1)));

      max_lambda_44 = shock_filter_control*fmax(fabs(-u2_B0(0,0,0) + a_B0(0,0,0)), fabs(-u2_B0(0,0,1) + a_B0(0,0,1)));

       beta_0 = ((1.0/4.0))*((CS_02*max_lambda_00 + CF_02)*(CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (CS_02*max_lambda_00 + CF_02))*(CS_01*max_lambda_00
            + CF_01);

       beta_1 = ((1.0/4.0))*((CS_01*max_lambda_00 + CF_01)*(CS_01*max_lambda_00 + CF_01)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - (CS_01*max_lambda_00 + CF_01))*(CS_00*max_lambda_00
            + CF_00);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_0 = (-(1.0/4.0)*(CS_00*max_lambda_00 + CF_00) + ((3.0/4.0))*(CS_01*max_lambda_00 + CF_01))*omega_1 +
            (((1.0/4.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(CS_02*max_lambda_00 + CF_02))*omega_0 + Recon_0;

       beta_0 = ((1.0/4.0))*((-CS_03*max_lambda_00 + CF_03)*(-CS_03*max_lambda_00 + CF_03)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (-CS_03*max_lambda_00 +
            CF_03))*(-CS_02*max_lambda_00 + CF_02);

       beta_1 = ((1.0/4.0))*((-CS_02*max_lambda_00 + CF_02)*(-CS_02*max_lambda_00 + CF_02)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - (-CS_02*max_lambda_00 +
            CF_02))*(-CS_01*max_lambda_00 + CF_01);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj0 = fmax(rj0, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_0 = (-(1.0/4.0)*(-CS_03*max_lambda_00 + CF_03) + ((3.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_0 +
            (((1.0/4.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/4.0))*(-CS_02*max_lambda_00 + CF_02))*omega_1 +
            Recon_0;

       beta_0 = ((1.0/4.0))*((CS_12*max_lambda_11 + CF_12)*(CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (CS_12*max_lambda_11 + CF_12))*(CS_11*max_lambda_11
            + CF_11);

       beta_1 = ((1.0/4.0))*((CS_11*max_lambda_11 + CF_11)*(CS_11*max_lambda_11 + CF_11)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - (CS_11*max_lambda_11 + CF_11))*(CS_10*max_lambda_11
            + CF_10);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_1 = (-(1.0/4.0)*(CS_10*max_lambda_11 + CF_10) + ((3.0/4.0))*(CS_11*max_lambda_11 + CF_11))*omega_1 +
            (((1.0/4.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(CS_12*max_lambda_11 + CF_12))*omega_0 + Recon_1;

       beta_0 = ((1.0/4.0))*((-CS_13*max_lambda_11 + CF_13)*(-CS_13*max_lambda_11 + CF_13)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (-CS_13*max_lambda_11 +
            CF_13))*(-CS_12*max_lambda_11 + CF_12);

       beta_1 = ((1.0/4.0))*((-CS_12*max_lambda_11 + CF_12)*(-CS_12*max_lambda_11 + CF_12)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - (-CS_12*max_lambda_11 +
            CF_12))*(-CS_11*max_lambda_11 + CF_11);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj1 = fmax(rj1, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_1 = (-(1.0/4.0)*(-CS_13*max_lambda_11 + CF_13) + ((3.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_0 +
            (((1.0/4.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/4.0))*(-CS_12*max_lambda_11 + CF_12))*omega_1 +
            Recon_1;

       beta_0 = ((1.0/4.0))*((CS_22*max_lambda_22 + CF_22)*(CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (CS_22*max_lambda_22 + CF_22))*(CS_21*max_lambda_22
            + CF_21);

       beta_1 = ((1.0/4.0))*((CS_21*max_lambda_22 + CF_21)*(CS_21*max_lambda_22 + CF_21)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - (CS_21*max_lambda_22 + CF_21))*(CS_20*max_lambda_22
            + CF_20);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_2 = (-(1.0/4.0)*(CS_20*max_lambda_22 + CF_20) + ((3.0/4.0))*(CS_21*max_lambda_22 + CF_21))*omega_1 +
            (((1.0/4.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(CS_22*max_lambda_22 + CF_22))*omega_0 + Recon_2;

       beta_0 = ((1.0/4.0))*((-CS_23*max_lambda_22 + CF_23)*(-CS_23*max_lambda_22 + CF_23)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (-CS_23*max_lambda_22 +
            CF_23))*(-CS_22*max_lambda_22 + CF_22);

       beta_1 = ((1.0/4.0))*((-CS_22*max_lambda_22 + CF_22)*(-CS_22*max_lambda_22 + CF_22)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - (-CS_22*max_lambda_22 +
            CF_22))*(-CS_21*max_lambda_22 + CF_21);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj2 = fmax(rj2, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_2 = (-(1.0/4.0)*(-CS_23*max_lambda_22 + CF_23) + ((3.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_0 +
            (((1.0/4.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/4.0))*(-CS_22*max_lambda_22 + CF_22))*omega_1 +
            Recon_2;

       beta_0 = ((1.0/4.0))*((CS_32*max_lambda_33 + CF_32)*(CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (CS_32*max_lambda_33 + CF_32))*(CS_31*max_lambda_33
            + CF_31);

       beta_1 = ((1.0/4.0))*((CS_31*max_lambda_33 + CF_31)*(CS_31*max_lambda_33 + CF_31)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - (CS_31*max_lambda_33 + CF_31))*(CS_30*max_lambda_33
            + CF_30);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_3 = (-(1.0/4.0)*(CS_30*max_lambda_33 + CF_30) + ((3.0/4.0))*(CS_31*max_lambda_33 + CF_31))*omega_1 +
            (((1.0/4.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(CS_32*max_lambda_33 + CF_32))*omega_0 + Recon_3;

       beta_0 = ((1.0/4.0))*((-CS_33*max_lambda_33 + CF_33)*(-CS_33*max_lambda_33 + CF_33)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (-CS_33*max_lambda_33 +
            CF_33))*(-CS_32*max_lambda_33 + CF_32);

       beta_1 = ((1.0/4.0))*((-CS_32*max_lambda_33 + CF_32)*(-CS_32*max_lambda_33 + CF_32)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - (-CS_32*max_lambda_33 +
            CF_32))*(-CS_31*max_lambda_33 + CF_31);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj3 = fmax(rj3, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_3 = (-(1.0/4.0)*(-CS_33*max_lambda_33 + CF_33) + ((3.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_0 +
            (((1.0/4.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/4.0))*(-CS_32*max_lambda_33 + CF_32))*omega_1 +
            Recon_3;

       beta_0 = ((1.0/4.0))*((CS_42*max_lambda_44 + CF_42)*(CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (CS_42*max_lambda_44 + CF_42))*(CS_41*max_lambda_44
            + CF_41);

       beta_1 = ((1.0/4.0))*((CS_41*max_lambda_44 + CF_41)*(CS_41*max_lambda_44 + CF_41)) +
            ((1.0/2.0))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - (CS_41*max_lambda_44 + CF_41))*(CS_40*max_lambda_44
            + CF_40);

       alpha_0 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = 0.333333333333333*fabs(-1.0 + 3*omega_1) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_0);

       Recon_4 = (-(1.0/4.0)*(CS_40*max_lambda_44 + CF_40) + ((3.0/4.0))*(CS_41*max_lambda_44 + CF_41))*omega_1 +
            (((1.0/4.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(CS_42*max_lambda_44 + CF_42))*omega_0 + Recon_4;

       beta_0 = ((1.0/4.0))*((-CS_43*max_lambda_44 + CF_43)*(-CS_43*max_lambda_44 + CF_43)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (-CS_43*max_lambda_44 +
            CF_43))*(-CS_42*max_lambda_44 + CF_42);

       beta_1 = ((1.0/4.0))*((-CS_42*max_lambda_44 + CF_42)*(-CS_42*max_lambda_44 + CF_42)) +
            ((1.0/2.0))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - (-CS_42*max_lambda_44 +
            CF_42))*(-CS_41*max_lambda_44 + CF_41);

       alpha_0 = 0.333333333333333 + ((1.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_0)*(1.0e-40 + beta_0));

       alpha_1 = 0.666666666666667 + ((2.0/3.0))*(fabs(-beta_1 + beta_0)*fabs(-beta_1 + beta_0))/((1.0e-40 +
            beta_1)*(1.0e-40 + beta_1));

      inv_alpha_sum = 1.0/((alpha_0 + alpha_1));

      omega_0 = alpha_0*inv_alpha_sum;

      omega_1 = alpha_1*inv_alpha_sum;

      rj4 = fmax(rj4, 0.333333333333333*fabs(-1.0 + 3*omega_0) + 0.333333333333333*fabs(-1.0 + ((3.0/2.0))*omega_1));

       Recon_4 = (-(1.0/4.0)*(-CS_43*max_lambda_44 + CF_43) + ((3.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_0 +
            (((1.0/4.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/4.0))*(-CS_42*max_lambda_44 + CF_42))*omega_1 +
            Recon_4;

      Recon_0 = (((1.0/12.0))*(-7*CF_01 - 7*CF_02 + CF_00 + CF_03) + Recon_0)*rj0;

      Recon_1 = (((1.0/12.0))*(-7*CF_11 - 7*CF_12 + CF_10 + CF_13) + Recon_1)*rj1;

      Recon_2 = (((1.0/12.0))*(-7*CF_21 - 7*CF_22 + CF_20 + CF_23) + Recon_2)*rj2;

      Recon_3 = (((1.0/12.0))*(-7*CF_31 - 7*CF_32 + CF_30 + CF_33) + Recon_3)*rj3;

      Recon_4 = (((1.0/12.0))*(-7*CF_41 - 7*CF_42 + CF_40 + CF_43) + Recon_4)*rj4;

       rho_RKold_B0(0,0,0) = 0.707106781186547*AVG_2_rho*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*Recon_4*inv_AVG_a + Recon_2;

       rhou0_RKold_B0(0,0,0) = AVG_2_u0*Recon_2 - AVG_2_rho*Recon_1 +
            0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_4*inv_AVG_a;

       rhou1_RKold_B0(0,0,0) = AVG_2_rho*Recon_0 + AVG_2_u1*Recon_2 +
            0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_3*inv_AVG_a +
            0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_4*inv_AVG_a;

       rhou2_RKold_B0(0,0,0) = AVG_2_u2*Recon_2 + 0.707106781186547*(-AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a +
            0.707106781186547*(AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a;

       rhoE_RKold_B0(0,0,0) = (((1.0/2.0))*(AVG_2_u0*AVG_2_u0) + ((1.0/2.0))*(AVG_2_u1*AVG_2_u1) +
            ((1.0/2.0))*(AVG_2_u2*AVG_2_u2))*Recon_2 + AVG_2_rho*AVG_2_u1*Recon_0 - AVG_2_rho*AVG_2_u0*Recon_1 +
            0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
            (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 + AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a +
            0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
            (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 - AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a;

   }

   else{

      rho_RKold_B0(0,0,0) = 0.0;

      rhou0_RKold_B0(0,0,0) = 0.0;

      rhou1_RKold_B0(0,0,0) = 0.0;

      rhou2_RKold_B0(0,0,0) = 0.0;

      rhoE_RKold_B0(0,0,0) = 0.0;

   }

}

 void opensbliblock00Kernel049(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
&Residual2_B0, const ACC<double> &Residual3_B0, const ACC<double> &Residual4_B0, const ACC<double> &kappa_B0, const
ACC<double> &rhoE_RKold_B0, const ACC<double> &rho_RKold_B0, const ACC<double> &rhou0_RKold_B0, const ACC<double>
&rhou1_RKold_B0, const ACC<double> &rhou2_RKold_B0, const ACC<double> &wk0_B0, const ACC<double> &wk1_B0, const
ACC<double> &wk2_B0, const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, ACC<double> &WENO_filter_B0, ACC<double>
&rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0, ACC<double> &rhou2_B0, const int *idx)
{
   double Grid_0 = 0.0;
   double Grid_1 = 0.0;
   double Grid_2 = 0.0;
   double Wall = 0.0;
   Grid_0 = idx[0];

   Grid_1 = idx[1];

   Grid_2 = idx[2];

   Wall = 1;

    WENO_filter_B0(0,0,0) = ((fmax(kappa_B0(0,0,1), fmax(kappa_B0(0,0,0), fmax(kappa_B0(0,2,0), fmax(kappa_B0(1,0,0),
      fmax(kappa_B0(0,-1,0), fmax(kappa_B0(-1,0,0), fmax(kappa_B0(0,1,0), fmax(kappa_B0(0,0,2), fmax(kappa_B0(2,0,0),
      kappa_B0(0,0,-1)))))))))) >= Ducros_select) ? (
   1
)
: (
   0.0
));

    rho_B0(0,0,0) = (-(-wk0_B0(-1,0,0) + wk0_B0(0,0,0))*inv_rfact0_block0 - (-rho_RKold_B0(0,0,-1) +
      rho_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual0_B0(0,-1,0) +
      Residual0_B0(0,0,0))*inv_rfact1_block0*Wall)*dt*WENO_filter_B0(0,0,0) + rho_B0(0,0,0);

    rhou0_B0(0,0,0) = (-(-wk1_B0(-1,0,0) + wk1_B0(0,0,0))*inv_rfact0_block0 - (-rhou0_RKold_B0(0,0,-1) +
      rhou0_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual1_B0(0,-1,0) +
      Residual1_B0(0,0,0))*inv_rfact1_block0*Wall)*dt*WENO_filter_B0(0,0,0) + rhou0_B0(0,0,0);

    rhou1_B0(0,0,0) = (-(-wk2_B0(-1,0,0) + wk2_B0(0,0,0))*inv_rfact0_block0 - (-rhou1_RKold_B0(0,0,-1) +
      rhou1_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual2_B0(0,-1,0) +
      Residual2_B0(0,0,0))*inv_rfact1_block0*Wall)*dt*WENO_filter_B0(0,0,0) + rhou1_B0(0,0,0);

    rhou2_B0(0,0,0) = (-(-wk3_B0(-1,0,0) + wk3_B0(0,0,0))*inv_rfact0_block0 - (-rhou2_RKold_B0(0,0,-1) +
      rhou2_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual3_B0(0,-1,0) +
      Residual3_B0(0,0,0))*inv_rfact1_block0*Wall)*dt*WENO_filter_B0(0,0,0) + rhou2_B0(0,0,0);

    rhoE_B0(0,0,0) = (-(-wk4_B0(-1,0,0) + wk4_B0(0,0,0))*inv_rfact0_block0 - (-rhoE_RKold_B0(0,0,-1) +
      rhoE_RKold_B0(0,0,0))*inv_rfact2_block0 - (-Residual4_B0(0,-1,0) +
      Residual4_B0(0,0,0))*inv_rfact1_block0*Wall)*dt*WENO_filter_B0(0,0,0) + rhoE_B0(0,0,0);

}

 void opensbliblock00Kernel042(const ACC<double> &mu_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, double *KE_B0, double *dilatation_dissipation_B0, double
*enstrophy_dissipation_B0, double *rhom_B0, ACC<double> &divV_B0)
{
   double d1_u0_dx = 0.0;
   double d1_u0_dy = 0.0;
   double d1_u0_dz = 0.0;
   double d1_u1_dx = 0.0;
   double d1_u1_dy = 0.0;
   double d1_u1_dz = 0.0;
   double d1_u2_dx = 0.0;
   double d1_u2_dy = 0.0;
   double d1_u2_dz = 0.0;
   double wx = 0.0;
   double wy = 0.0;
   double wz = 0.0;
    d1_u2_dy = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
      ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

    d1_u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

   wx = -d1_u1_dz + d1_u2_dy;

    d1_u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

    d1_u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

   wy = -d1_u2_dx + d1_u0_dz;

    d1_u0_dy = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
      ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

    d1_u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

   wz = -d1_u0_dy + d1_u1_dx;

    d1_u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

    d1_u1_dy = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
      ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

    d1_u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

   divV_B0(0,0,0) = d1_u0_dx + d1_u1_dy + d1_u2_dz;

   *rhom_B0 = rho_B0(0,0,0) + *rhom_B0;

    *KE_B0 = 0.5*((u0_B0(0,0,0)*u0_B0(0,0,0)) + (u1_B0(0,0,0)*u1_B0(0,0,0)) + (u2_B0(0,0,0)*u2_B0(0,0,0)))*rho_B0(0,0,0)
      + *KE_B0;

   *dilatation_dissipation_B0 = ((4.0/3.0))*(divV_B0(0,0,0)*divV_B0(0,0,0))*mu_B0(0,0,0) + *dilatation_dissipation_B0;

   *enstrophy_dissipation_B0 = ((wx*wx) + (wy*wy) + (wz*wz))*mu_B0(0,0,0) + *enstrophy_dissipation_B0;

}

#endif
