#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
 void opensbliblock00Kernel043(ACC<double> &rhoE_B0, ACC<double> &rho_B0, ACC<double> &rhou0_B0, ACC<double> &rhou1_B0,
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

void opensbliblock00Kernel003(const ACC<double> &rho_B0, const ACC<double> &rhou0_B0, ACC<double> &u0_B0)
{
   u0_B0(0,0,0) = rhou0_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel011(const ACC<double> &rho_B0, const ACC<double> &rhou2_B0, ACC<double> &u2_B0)
{
   u2_B0(0,0,0) = rhou2_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel012(const ACC<double> &rho_B0, const ACC<double> &rhou1_B0, ACC<double> &u1_B0)
{
   u1_B0(0,0,0) = rhou1_B0(0,0,0)/rho_B0(0,0,0);

}

 void opensbliblock00Kernel005(const ACC<double> &rhoE_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &p_B0)
{
    p_B0(0,0,0) = (-1 + gama)*(-(1.0/2.0)*(u0_B0(0,0,0)*u0_B0(0,0,0))*rho_B0(0,0,0) -
      (1.0/2.0)*(u1_B0(0,0,0)*u1_B0(0,0,0))*rho_B0(0,0,0) - (1.0/2.0)*(u2_B0(0,0,0)*u2_B0(0,0,0))*rho_B0(0,0,0) +
      rhoE_B0(0,0,0));

}

void opensbliblock00Kernel006(const ACC<double> &p_B0, const ACC<double> &rho_B0, ACC<double> &a_B0)
{
   a_B0(0,0,0) = sqrt((gama*p_B0(0,0,0)/rho_B0(0,0,0)));

}

void opensbliblock00Kernel015(const ACC<double> &p_B0, const ACC<double> &rho_B0, ACC<double> &T_B0)
{
   T_B0(0,0,0) = (Minf*Minf)*gama*p_B0(0,0,0)/rho_B0(0,0,0);

}

void opensbliblock00Kernel014(const ACC<double> &T_B0, ACC<double> &mu_B0)
{
   mu_B0(0,0,0) = 1.4042*T_B0(0,0,0)*sqrt(T_B0(0,0,0))/(0.40417 + T_B0(0,0,0));

}

 void opensbliblock00Kernel000(const ACC<double> &a_B0, const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const
ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, const
ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &wk0_B0, ACC<double> &wk1_B0,
ACC<double> &wk2_B0, ACC<double> &wk3_B0, ACC<double> &wk4_B0)
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
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double delta_0 = 0.0;
   double delta_1 = 0.0;
   double delta_2 = 0.0;
   double delta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double inv_beta_0 = 0.0;
   double inv_beta_1 = 0.0;
   double inv_beta_2 = 0.0;
   double inv_beta_3 = 0.0;
   double inv_omega_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
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

    CF_00 = p_B0(-2,0,0)*AVG_0_0_LEV_01 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_00 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_04 +
      u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_04 + u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_01 +
      u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_02 + u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_03;

   CF_10 = rhou0_B0(-2,0,0)*AVG_0_0_LEV_10 + u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_13;

   CF_20 = rhou0_B0(-2,0,0)*AVG_0_0_LEV_20 + u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_22;

    CF_30 = p_B0(-2,0,0)*AVG_0_0_LEV_31 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_30 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_34 +
      u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_34 + u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_31 +
      u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_32 + u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_33;

    CF_40 = p_B0(-2,0,0)*AVG_0_0_LEV_41 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_40 + p_B0(-2,0,0)*u0_B0(-2,0,0)*AVG_0_0_LEV_44 +
      u0_B0(-2,0,0)*rhoE_B0(-2,0,0)*AVG_0_0_LEV_44 + u0_B0(-2,0,0)*rhou0_B0(-2,0,0)*AVG_0_0_LEV_41 +
      u0_B0(-2,0,0)*rhou1_B0(-2,0,0)*AVG_0_0_LEV_42 + u0_B0(-2,0,0)*rhou2_B0(-2,0,0)*AVG_0_0_LEV_43;

    CS_00 = rho_B0(-2,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(-2,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_03;

   CS_10 = rho_B0(-2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_13;

   CS_20 = rho_B0(-2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-2,0,0)*AVG_0_0_LEV_22;

    CS_30 = rho_B0(-2,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(-2,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_33;

    CS_40 = rho_B0(-2,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-2,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-2,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(-2,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-2,0,0)*AVG_0_0_LEV_43;

    CF_01 = p_B0(-1,0,0)*AVG_0_0_LEV_01 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_00 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_04 +
      u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 +
      u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

   CF_11 = rhou0_B0(-1,0,0)*AVG_0_0_LEV_10 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_13;

   CF_21 = rhou0_B0(-1,0,0)*AVG_0_0_LEV_20 + u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_22;

    CF_31 = p_B0(-1,0,0)*AVG_0_0_LEV_31 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_30 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_34 +
      u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 +
      u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

    CF_41 = p_B0(-1,0,0)*AVG_0_0_LEV_41 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_40 + p_B0(-1,0,0)*u0_B0(-1,0,0)*AVG_0_0_LEV_44 +
      u0_B0(-1,0,0)*rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 + u0_B0(-1,0,0)*rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 +
      u0_B0(-1,0,0)*rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 + u0_B0(-1,0,0)*rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

    CS_01 = rho_B0(-1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(-1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_03;

   CS_11 = rho_B0(-1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_13;

   CS_21 = rho_B0(-1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(-1,0,0)*AVG_0_0_LEV_22;

    CS_31 = rho_B0(-1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(-1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_33;

    CS_41 = rho_B0(-1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(-1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(-1,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(-1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(-1,0,0)*AVG_0_0_LEV_43;

    CF_02 = p_B0(0,0,0)*AVG_0_0_LEV_01 + rhou0_B0(0,0,0)*AVG_0_0_LEV_00 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_04 +
      u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_04 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_01 +
      u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_02 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

   CF_12 = rhou0_B0(0,0,0)*AVG_0_0_LEV_10 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_13;

   CF_22 = rhou0_B0(0,0,0)*AVG_0_0_LEV_20 + u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_22;

    CF_32 = p_B0(0,0,0)*AVG_0_0_LEV_31 + rhou0_B0(0,0,0)*AVG_0_0_LEV_30 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_34 +
      u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_34 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_31 +
      u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_32 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

    CF_42 = p_B0(0,0,0)*AVG_0_0_LEV_41 + rhou0_B0(0,0,0)*AVG_0_0_LEV_40 + p_B0(0,0,0)*u0_B0(0,0,0)*AVG_0_0_LEV_44 +
      u0_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_0_0_LEV_44 + u0_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_0_0_LEV_41 +
      u0_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_0_0_LEV_42 + u0_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

    CS_02 = rho_B0(0,0,0)*AVG_0_0_LEV_00 + rhoE_B0(0,0,0)*AVG_0_0_LEV_04 + rhou0_B0(0,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(0,0,0)*AVG_0_0_LEV_02 + rhou2_B0(0,0,0)*AVG_0_0_LEV_03;

   CS_12 = rho_B0(0,0,0)*AVG_0_0_LEV_10 + rhou2_B0(0,0,0)*AVG_0_0_LEV_13;

   CS_22 = rho_B0(0,0,0)*AVG_0_0_LEV_20 + rhou1_B0(0,0,0)*AVG_0_0_LEV_22;

    CS_32 = rho_B0(0,0,0)*AVG_0_0_LEV_30 + rhoE_B0(0,0,0)*AVG_0_0_LEV_34 + rhou0_B0(0,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(0,0,0)*AVG_0_0_LEV_32 + rhou2_B0(0,0,0)*AVG_0_0_LEV_33;

    CS_42 = rho_B0(0,0,0)*AVG_0_0_LEV_40 + rhoE_B0(0,0,0)*AVG_0_0_LEV_44 + rhou0_B0(0,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(0,0,0)*AVG_0_0_LEV_42 + rhou2_B0(0,0,0)*AVG_0_0_LEV_43;

    CF_03 = p_B0(1,0,0)*AVG_0_0_LEV_01 + rhou0_B0(1,0,0)*AVG_0_0_LEV_00 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_04 +
      u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_04 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_01 +
      u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_02 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

   CF_13 = rhou0_B0(1,0,0)*AVG_0_0_LEV_10 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_13;

   CF_23 = rhou0_B0(1,0,0)*AVG_0_0_LEV_20 + u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_22;

    CF_33 = p_B0(1,0,0)*AVG_0_0_LEV_31 + rhou0_B0(1,0,0)*AVG_0_0_LEV_30 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_34 +
      u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_34 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_31 +
      u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_32 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

    CF_43 = p_B0(1,0,0)*AVG_0_0_LEV_41 + rhou0_B0(1,0,0)*AVG_0_0_LEV_40 + p_B0(1,0,0)*u0_B0(1,0,0)*AVG_0_0_LEV_44 +
      u0_B0(1,0,0)*rhoE_B0(1,0,0)*AVG_0_0_LEV_44 + u0_B0(1,0,0)*rhou0_B0(1,0,0)*AVG_0_0_LEV_41 +
      u0_B0(1,0,0)*rhou1_B0(1,0,0)*AVG_0_0_LEV_42 + u0_B0(1,0,0)*rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

    CS_03 = rho_B0(1,0,0)*AVG_0_0_LEV_00 + rhoE_B0(1,0,0)*AVG_0_0_LEV_04 + rhou0_B0(1,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(1,0,0)*AVG_0_0_LEV_02 + rhou2_B0(1,0,0)*AVG_0_0_LEV_03;

   CS_13 = rho_B0(1,0,0)*AVG_0_0_LEV_10 + rhou2_B0(1,0,0)*AVG_0_0_LEV_13;

   CS_23 = rho_B0(1,0,0)*AVG_0_0_LEV_20 + rhou1_B0(1,0,0)*AVG_0_0_LEV_22;

    CS_33 = rho_B0(1,0,0)*AVG_0_0_LEV_30 + rhoE_B0(1,0,0)*AVG_0_0_LEV_34 + rhou0_B0(1,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(1,0,0)*AVG_0_0_LEV_32 + rhou2_B0(1,0,0)*AVG_0_0_LEV_33;

    CS_43 = rho_B0(1,0,0)*AVG_0_0_LEV_40 + rhoE_B0(1,0,0)*AVG_0_0_LEV_44 + rhou0_B0(1,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(1,0,0)*AVG_0_0_LEV_42 + rhou2_B0(1,0,0)*AVG_0_0_LEV_43;

    CF_04 = p_B0(2,0,0)*AVG_0_0_LEV_01 + rhou0_B0(2,0,0)*AVG_0_0_LEV_00 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_04 +
      u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_04 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_01 +
      u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_02 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

   CF_14 = rhou0_B0(2,0,0)*AVG_0_0_LEV_10 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_13;

   CF_24 = rhou0_B0(2,0,0)*AVG_0_0_LEV_20 + u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_22;

    CF_34 = p_B0(2,0,0)*AVG_0_0_LEV_31 + rhou0_B0(2,0,0)*AVG_0_0_LEV_30 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_34 +
      u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_34 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_31 +
      u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_32 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

    CF_44 = p_B0(2,0,0)*AVG_0_0_LEV_41 + rhou0_B0(2,0,0)*AVG_0_0_LEV_40 + p_B0(2,0,0)*u0_B0(2,0,0)*AVG_0_0_LEV_44 +
      u0_B0(2,0,0)*rhoE_B0(2,0,0)*AVG_0_0_LEV_44 + u0_B0(2,0,0)*rhou0_B0(2,0,0)*AVG_0_0_LEV_41 +
      u0_B0(2,0,0)*rhou1_B0(2,0,0)*AVG_0_0_LEV_42 + u0_B0(2,0,0)*rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

    CS_04 = rho_B0(2,0,0)*AVG_0_0_LEV_00 + rhoE_B0(2,0,0)*AVG_0_0_LEV_04 + rhou0_B0(2,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(2,0,0)*AVG_0_0_LEV_02 + rhou2_B0(2,0,0)*AVG_0_0_LEV_03;

   CS_14 = rho_B0(2,0,0)*AVG_0_0_LEV_10 + rhou2_B0(2,0,0)*AVG_0_0_LEV_13;

   CS_24 = rho_B0(2,0,0)*AVG_0_0_LEV_20 + rhou1_B0(2,0,0)*AVG_0_0_LEV_22;

    CS_34 = rho_B0(2,0,0)*AVG_0_0_LEV_30 + rhoE_B0(2,0,0)*AVG_0_0_LEV_34 + rhou0_B0(2,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(2,0,0)*AVG_0_0_LEV_32 + rhou2_B0(2,0,0)*AVG_0_0_LEV_33;

    CS_44 = rho_B0(2,0,0)*AVG_0_0_LEV_40 + rhoE_B0(2,0,0)*AVG_0_0_LEV_44 + rhou0_B0(2,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(2,0,0)*AVG_0_0_LEV_42 + rhou2_B0(2,0,0)*AVG_0_0_LEV_43;

    CF_05 = p_B0(3,0,0)*AVG_0_0_LEV_01 + rhou0_B0(3,0,0)*AVG_0_0_LEV_00 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_04 +
      u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_04 + u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_01 +
      u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_02 + u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_03;

   CF_15 = rhou0_B0(3,0,0)*AVG_0_0_LEV_10 + u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_13;

   CF_25 = rhou0_B0(3,0,0)*AVG_0_0_LEV_20 + u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_22;

    CF_35 = p_B0(3,0,0)*AVG_0_0_LEV_31 + rhou0_B0(3,0,0)*AVG_0_0_LEV_30 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_34 +
      u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_34 + u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_31 +
      u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_32 + u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_33;

    CF_45 = p_B0(3,0,0)*AVG_0_0_LEV_41 + rhou0_B0(3,0,0)*AVG_0_0_LEV_40 + p_B0(3,0,0)*u0_B0(3,0,0)*AVG_0_0_LEV_44 +
      u0_B0(3,0,0)*rhoE_B0(3,0,0)*AVG_0_0_LEV_44 + u0_B0(3,0,0)*rhou0_B0(3,0,0)*AVG_0_0_LEV_41 +
      u0_B0(3,0,0)*rhou1_B0(3,0,0)*AVG_0_0_LEV_42 + u0_B0(3,0,0)*rhou2_B0(3,0,0)*AVG_0_0_LEV_43;

    CS_05 = rho_B0(3,0,0)*AVG_0_0_LEV_00 + rhoE_B0(3,0,0)*AVG_0_0_LEV_04 + rhou0_B0(3,0,0)*AVG_0_0_LEV_01 +
      rhou1_B0(3,0,0)*AVG_0_0_LEV_02 + rhou2_B0(3,0,0)*AVG_0_0_LEV_03;

   CS_15 = rho_B0(3,0,0)*AVG_0_0_LEV_10 + rhou2_B0(3,0,0)*AVG_0_0_LEV_13;

   CS_25 = rho_B0(3,0,0)*AVG_0_0_LEV_20 + rhou1_B0(3,0,0)*AVG_0_0_LEV_22;

    CS_35 = rho_B0(3,0,0)*AVG_0_0_LEV_30 + rhoE_B0(3,0,0)*AVG_0_0_LEV_34 + rhou0_B0(3,0,0)*AVG_0_0_LEV_31 +
      rhou1_B0(3,0,0)*AVG_0_0_LEV_32 + rhou2_B0(3,0,0)*AVG_0_0_LEV_33;

    CS_45 = rho_B0(3,0,0)*AVG_0_0_LEV_40 + rhoE_B0(3,0,0)*AVG_0_0_LEV_44 + rhou0_B0(3,0,0)*AVG_0_0_LEV_41 +
      rhou1_B0(3,0,0)*AVG_0_0_LEV_42 + rhou2_B0(3,0,0)*AVG_0_0_LEV_43;

   max_lambda_00 = shock_filter_control*fmax(fabs(u0_B0(0,0,0)), fabs(u0_B0(1,0,0)));

   max_lambda_11 = max_lambda_00;

   max_lambda_22 = max_lambda_00;

   max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(1,0,0) + u0_B0(1,0,0)), fabs(a_B0(0,0,0) + u0_B0(0,0,0)));

   max_lambda_44 = shock_filter_control*fmax(fabs(-u0_B0(1,0,0) + a_B0(1,0,0)), fabs(-u0_B0(0,0,0) + a_B0(0,0,0)));

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 +
      CF_03))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 + CF_03))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03) -
      (CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 +
      CF_03) - (CS_02*max_lambda_00 + CF_02)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 + CF_03) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 +
      CF_03) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(CS_02*max_lambda_00 +
      CF_02) + ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 + CF_01) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 +
      CF_01) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) +
      ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01))*(((1.0/2.0))*(CS_00*max_lambda_00 +
      CF_00) + ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01)));

    beta_3 = -(781.0/480.0)*(CS_04*max_lambda_00 + CF_04) - (781.0/1440.0)*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/36.0))*((9*(CS_03*max_lambda_00 + CF_03) - (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) -
      (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05)*(9*(CS_03*max_lambda_00 + CF_03) -
      (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) - (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05))
      + ((13.0/12.0))*((2*(CS_04*max_lambda_00 + CF_04) - (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) -
      (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02)*(2*(CS_04*max_lambda_00 + CF_04) -
      (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) - (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02))
      + ((781.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((781.0/1440.0))*(CS_05*max_lambda_00 + CF_05);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(CS_03*max_lambda_00 + CF_03))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/8.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/24.0))*(CS_05*max_lambda_00 + CF_05) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_00*max_lambda_00 + CF_00) +
      ((11.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 + CF_04))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04) -
      (-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
      CF_04) - (-CS_03*max_lambda_00 + CF_03)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00 + CF_02) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00
      + CF_02) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (-CS_02*max_lambda_00 +
      CF_02))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) -
      (-CS_02*max_lambda_00 + CF_02)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00 + CF_04) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00
      + CF_04) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 +
      CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - (-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) -
      (-CS_04*max_lambda_00 + CF_04)));

    beta_3 = ((1.0/36.0))*((-(-CS_00*max_lambda_00 + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) +
      ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(-(-CS_00*max_lambda_00
      + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) + ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) +
      ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((2*(-CS_01*max_lambda_00 + CF_01) -
      (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 +
      CF_03)*(2*(-CS_01*max_lambda_00 + CF_01) - (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 + CF_03)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 + CF_01)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(-CS_02*max_lambda_00 + CF_02))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((1.0/24.0))*(-CS_00*max_lambda_00 + CF_00) + ((13.0/24.0))*(-CS_02*max_lambda_00 + CF_02))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_05*max_lambda_00 + CF_05) +
      ((11.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 +
      CF_13))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 + CF_13))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13) -
      (CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 +
      CF_13) - (CS_12*max_lambda_11 + CF_12)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 + CF_13) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 +
      CF_13) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(CS_12*max_lambda_11 +
      CF_12) + ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 + CF_11) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 +
      CF_11) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) +
      ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11))*(((1.0/2.0))*(CS_10*max_lambda_11 +
      CF_10) + ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11)));

    beta_3 = -(781.0/480.0)*(CS_14*max_lambda_11 + CF_14) - (781.0/1440.0)*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/36.0))*((9*(CS_13*max_lambda_11 + CF_13) - (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) -
      (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15)*(9*(CS_13*max_lambda_11 + CF_13) -
      (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) - (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15))
      + ((13.0/12.0))*((2*(CS_14*max_lambda_11 + CF_14) - (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) -
      (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12)*(2*(CS_14*max_lambda_11 + CF_14) -
      (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) - (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12))
      + ((781.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((781.0/1440.0))*(CS_15*max_lambda_11 + CF_15);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(CS_13*max_lambda_11 + CF_13))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/8.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/24.0))*(CS_15*max_lambda_11 + CF_15) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_10*max_lambda_11 + CF_10) +
      ((11.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 + CF_14))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14) -
      (-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
      CF_14) - (-CS_13*max_lambda_11 + CF_13)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11 + CF_12) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11
      + CF_12) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (-CS_12*max_lambda_11 +
      CF_12))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) -
      (-CS_12*max_lambda_11 + CF_12)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11 + CF_14) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11
      + CF_14) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 +
      CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - (-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) -
      (-CS_14*max_lambda_11 + CF_14)));

    beta_3 = ((1.0/36.0))*((-(-CS_10*max_lambda_11 + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) +
      ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(-(-CS_10*max_lambda_11
      + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) + ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) +
      ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((2*(-CS_11*max_lambda_11 + CF_11) -
      (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 +
      CF_13)*(2*(-CS_11*max_lambda_11 + CF_11) - (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 + CF_13)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 + CF_11)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(-CS_12*max_lambda_11 + CF_12))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((1.0/24.0))*(-CS_10*max_lambda_11 + CF_10) + ((13.0/24.0))*(-CS_12*max_lambda_11 + CF_12))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_15*max_lambda_11 + CF_15) +
      ((11.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 +
      CF_23))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 + CF_23))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23) -
      (CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 +
      CF_23) - (CS_22*max_lambda_22 + CF_22)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 + CF_23) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 +
      CF_23) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(CS_22*max_lambda_22 +
      CF_22) + ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 + CF_21) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 +
      CF_21) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) +
      ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21))*(((1.0/2.0))*(CS_20*max_lambda_22 +
      CF_20) + ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21)));

    beta_3 = -(781.0/480.0)*(CS_24*max_lambda_22 + CF_24) - (781.0/1440.0)*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/36.0))*((9*(CS_23*max_lambda_22 + CF_23) - (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) -
      (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25)*(9*(CS_23*max_lambda_22 + CF_23) -
      (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) - (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25))
      + ((13.0/12.0))*((2*(CS_24*max_lambda_22 + CF_24) - (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) -
      (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22)*(2*(CS_24*max_lambda_22 + CF_24) -
      (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) - (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22))
      + ((781.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((781.0/1440.0))*(CS_25*max_lambda_22 + CF_25);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(CS_23*max_lambda_22 + CF_23))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/8.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/24.0))*(CS_25*max_lambda_22 + CF_25) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_20*max_lambda_22 + CF_20) +
      ((11.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 + CF_24))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24) -
      (-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
      CF_24) - (-CS_23*max_lambda_22 + CF_23)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22 + CF_22) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22
      + CF_22) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (-CS_22*max_lambda_22 +
      CF_22))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) -
      (-CS_22*max_lambda_22 + CF_22)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22 + CF_24) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22
      + CF_24) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 +
      CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - (-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) -
      (-CS_24*max_lambda_22 + CF_24)));

    beta_3 = ((1.0/36.0))*((-(-CS_20*max_lambda_22 + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) +
      ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(-(-CS_20*max_lambda_22
      + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) + ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) +
      ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((2*(-CS_21*max_lambda_22 + CF_21) -
      (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 +
      CF_23)*(2*(-CS_21*max_lambda_22 + CF_21) - (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 + CF_23)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 + CF_21)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(-CS_22*max_lambda_22 + CF_22))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((1.0/24.0))*(-CS_20*max_lambda_22 + CF_20) + ((13.0/24.0))*(-CS_22*max_lambda_22 + CF_22))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_25*max_lambda_22 + CF_25) +
      ((11.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 +
      CF_33))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 + CF_33))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33) -
      (CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 +
      CF_33) - (CS_32*max_lambda_33 + CF_32)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 + CF_33) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 +
      CF_33) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(CS_32*max_lambda_33 +
      CF_32) + ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 + CF_31) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 +
      CF_31) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) +
      ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31))*(((1.0/2.0))*(CS_30*max_lambda_33 +
      CF_30) + ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31)));

    beta_3 = -(781.0/480.0)*(CS_34*max_lambda_33 + CF_34) - (781.0/1440.0)*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/36.0))*((9*(CS_33*max_lambda_33 + CF_33) - (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) -
      (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35)*(9*(CS_33*max_lambda_33 + CF_33) -
      (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) - (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35))
      + ((13.0/12.0))*((2*(CS_34*max_lambda_33 + CF_34) - (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) -
      (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32)*(2*(CS_34*max_lambda_33 + CF_34) -
      (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) - (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32))
      + ((781.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((781.0/1440.0))*(CS_35*max_lambda_33 + CF_35);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(CS_33*max_lambda_33 + CF_33))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/8.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/24.0))*(CS_35*max_lambda_33 + CF_35) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_30*max_lambda_33 + CF_30) +
      ((11.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 + CF_34))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34) -
      (-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
      CF_34) - (-CS_33*max_lambda_33 + CF_33)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33 + CF_32) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33
      + CF_32) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (-CS_32*max_lambda_33 +
      CF_32))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) -
      (-CS_32*max_lambda_33 + CF_32)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33 + CF_34) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33
      + CF_34) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 +
      CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - (-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) -
      (-CS_34*max_lambda_33 + CF_34)));

    beta_3 = ((1.0/36.0))*((-(-CS_30*max_lambda_33 + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) +
      ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(-(-CS_30*max_lambda_33
      + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) + ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) +
      ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((2*(-CS_31*max_lambda_33 + CF_31) -
      (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 +
      CF_33)*(2*(-CS_31*max_lambda_33 + CF_31) - (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 + CF_33)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 + CF_31)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(-CS_32*max_lambda_33 + CF_32))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((1.0/24.0))*(-CS_30*max_lambda_33 + CF_30) + ((13.0/24.0))*(-CS_32*max_lambda_33 + CF_32))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_35*max_lambda_33 + CF_35) +
      ((11.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 +
      CF_43))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 + CF_43))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43) -
      (CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 +
      CF_43) - (CS_42*max_lambda_44 + CF_42)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 + CF_43) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 +
      CF_43) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(CS_42*max_lambda_44 +
      CF_42) + ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 + CF_41) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 +
      CF_41) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) +
      ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41))*(((1.0/2.0))*(CS_40*max_lambda_44 +
      CF_40) + ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41)));

    beta_3 = -(781.0/480.0)*(CS_44*max_lambda_44 + CF_44) - (781.0/1440.0)*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/36.0))*((9*(CS_43*max_lambda_44 + CF_43) - (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) -
      (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45)*(9*(CS_43*max_lambda_44 + CF_43) -
      (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) - (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45))
      + ((13.0/12.0))*((2*(CS_44*max_lambda_44 + CF_44) - (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) -
      (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42)*(2*(CS_44*max_lambda_44 + CF_44) -
      (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) - (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42))
      + ((781.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((781.0/1440.0))*(CS_45*max_lambda_44 + CF_45);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(CS_43*max_lambda_44 + CF_43))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/8.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/24.0))*(CS_45*max_lambda_44 + CF_45) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_40*max_lambda_44 + CF_40) +
      ((11.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_0*inv_omega_sum + Recon_4;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 + CF_44))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44) -
      (-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
      CF_44) - (-CS_43*max_lambda_44 + CF_43)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44 + CF_42) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44
      + CF_42) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (-CS_42*max_lambda_44 +
      CF_42))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) -
      (-CS_42*max_lambda_44 + CF_42)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44 + CF_44) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44
      + CF_44) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 +
      CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - (-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) -
      (-CS_44*max_lambda_44 + CF_44)));

    beta_3 = ((1.0/36.0))*((-(-CS_40*max_lambda_44 + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) +
      ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(-(-CS_40*max_lambda_44
      + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) + ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) +
      ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((2*(-CS_41*max_lambda_44 + CF_41) -
      (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 +
      CF_43)*(2*(-CS_41*max_lambda_44 + CF_41) - (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 + CF_43)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 + CF_41)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(-CS_42*max_lambda_44 + CF_42))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((1.0/24.0))*(-CS_40*max_lambda_44 + CF_40) + ((13.0/24.0))*(-CS_42*max_lambda_44 + CF_42))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_45*max_lambda_44 + CF_45) +
      ((11.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_0*inv_omega_sum + Recon_4;

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

 void opensbliblock00Kernel001(const ACC<double> &a_B0, const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const
ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, const
ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &wk5_B0, ACC<double> &wk6_B0,
ACC<double> &wk7_B0, ACC<double> &wk8_B0, ACC<double> &wk9_B0)
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
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double delta_0 = 0.0;
   double delta_1 = 0.0;
   double delta_2 = 0.0;
   double delta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double inv_beta_0 = 0.0;
   double inv_beta_1 = 0.0;
   double inv_beta_2 = 0.0;
   double inv_beta_3 = 0.0;
   double inv_omega_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
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

   CF_00 = rhou1_B0(0,-2,0)*AVG_1_1_LEV_00 + u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_03;

    CF_10 = p_B0(0,-2,0)*AVG_1_1_LEV_12 + rhou1_B0(0,-2,0)*AVG_1_1_LEV_10 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_14 +
      u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_14 + u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_11 +
      u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_12 + u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_13;

   CF_20 = rhou1_B0(0,-2,0)*AVG_1_1_LEV_20 + u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_21;

    CF_30 = p_B0(0,-2,0)*AVG_1_1_LEV_32 + rhou1_B0(0,-2,0)*AVG_1_1_LEV_30 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_34 +
      u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_34 + u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_31 +
      u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_32 + u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_33;

    CF_40 = p_B0(0,-2,0)*AVG_1_1_LEV_42 + rhou1_B0(0,-2,0)*AVG_1_1_LEV_40 + p_B0(0,-2,0)*u1_B0(0,-2,0)*AVG_1_1_LEV_44 +
      u1_B0(0,-2,0)*rhoE_B0(0,-2,0)*AVG_1_1_LEV_44 + u1_B0(0,-2,0)*rhou0_B0(0,-2,0)*AVG_1_1_LEV_41 +
      u1_B0(0,-2,0)*rhou1_B0(0,-2,0)*AVG_1_1_LEV_42 + u1_B0(0,-2,0)*rhou2_B0(0,-2,0)*AVG_1_1_LEV_43;

   CS_00 = rho_B0(0,-2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_03;

    CS_10 = rho_B0(0,-2,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,-2,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_13;

   CS_20 = rho_B0(0,-2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_21;

    CS_30 = rho_B0(0,-2,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,-2,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_33;

    CS_40 = rho_B0(0,-2,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-2,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-2,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,-2,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-2,0)*AVG_1_1_LEV_43;

   CF_01 = rhou1_B0(0,-1,0)*AVG_1_1_LEV_00 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_03;

    CF_11 = p_B0(0,-1,0)*AVG_1_1_LEV_12 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_10 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_14 +
      u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 +
      u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

   CF_21 = rhou1_B0(0,-1,0)*AVG_1_1_LEV_20 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_21;

    CF_31 = p_B0(0,-1,0)*AVG_1_1_LEV_32 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_30 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_34 +
      u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 +
      u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

    CF_41 = p_B0(0,-1,0)*AVG_1_1_LEV_42 + rhou1_B0(0,-1,0)*AVG_1_1_LEV_40 + p_B0(0,-1,0)*u1_B0(0,-1,0)*AVG_1_1_LEV_44 +
      u1_B0(0,-1,0)*rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 + u1_B0(0,-1,0)*rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 +
      u1_B0(0,-1,0)*rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 + u1_B0(0,-1,0)*rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

   CS_01 = rho_B0(0,-1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_03;

    CS_11 = rho_B0(0,-1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,-1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_13;

   CS_21 = rho_B0(0,-1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_21;

    CS_31 = rho_B0(0,-1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,-1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_33;

    CS_41 = rho_B0(0,-1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,-1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,-1,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,-1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,-1,0)*AVG_1_1_LEV_43;

   CF_02 = rhou1_B0(0,0,0)*AVG_1_1_LEV_00 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_03;

    CF_12 = p_B0(0,0,0)*AVG_1_1_LEV_12 + rhou1_B0(0,0,0)*AVG_1_1_LEV_10 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_14 +
      u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_14 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_11 +
      u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_12 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

   CF_22 = rhou1_B0(0,0,0)*AVG_1_1_LEV_20 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_21;

    CF_32 = p_B0(0,0,0)*AVG_1_1_LEV_32 + rhou1_B0(0,0,0)*AVG_1_1_LEV_30 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_34 +
      u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_34 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_31 +
      u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_32 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

    CF_42 = p_B0(0,0,0)*AVG_1_1_LEV_42 + rhou1_B0(0,0,0)*AVG_1_1_LEV_40 + p_B0(0,0,0)*u1_B0(0,0,0)*AVG_1_1_LEV_44 +
      u1_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_1_1_LEV_44 + u1_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_1_1_LEV_41 +
      u1_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_1_1_LEV_42 + u1_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

   CS_02 = rho_B0(0,0,0)*AVG_1_1_LEV_00 + rhou2_B0(0,0,0)*AVG_1_1_LEV_03;

    CS_12 = rho_B0(0,0,0)*AVG_1_1_LEV_10 + rhoE_B0(0,0,0)*AVG_1_1_LEV_14 + rhou0_B0(0,0,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,0,0)*AVG_1_1_LEV_12 + rhou2_B0(0,0,0)*AVG_1_1_LEV_13;

   CS_22 = rho_B0(0,0,0)*AVG_1_1_LEV_20 + rhou0_B0(0,0,0)*AVG_1_1_LEV_21;

    CS_32 = rho_B0(0,0,0)*AVG_1_1_LEV_30 + rhoE_B0(0,0,0)*AVG_1_1_LEV_34 + rhou0_B0(0,0,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,0,0)*AVG_1_1_LEV_32 + rhou2_B0(0,0,0)*AVG_1_1_LEV_33;

    CS_42 = rho_B0(0,0,0)*AVG_1_1_LEV_40 + rhoE_B0(0,0,0)*AVG_1_1_LEV_44 + rhou0_B0(0,0,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,0,0)*AVG_1_1_LEV_42 + rhou2_B0(0,0,0)*AVG_1_1_LEV_43;

   CF_03 = rhou1_B0(0,1,0)*AVG_1_1_LEV_00 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_03;

    CF_13 = p_B0(0,1,0)*AVG_1_1_LEV_12 + rhou1_B0(0,1,0)*AVG_1_1_LEV_10 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_14 +
      u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_14 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_11 +
      u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_12 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

   CF_23 = rhou1_B0(0,1,0)*AVG_1_1_LEV_20 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_21;

    CF_33 = p_B0(0,1,0)*AVG_1_1_LEV_32 + rhou1_B0(0,1,0)*AVG_1_1_LEV_30 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_34 +
      u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_34 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_31 +
      u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_32 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

    CF_43 = p_B0(0,1,0)*AVG_1_1_LEV_42 + rhou1_B0(0,1,0)*AVG_1_1_LEV_40 + p_B0(0,1,0)*u1_B0(0,1,0)*AVG_1_1_LEV_44 +
      u1_B0(0,1,0)*rhoE_B0(0,1,0)*AVG_1_1_LEV_44 + u1_B0(0,1,0)*rhou0_B0(0,1,0)*AVG_1_1_LEV_41 +
      u1_B0(0,1,0)*rhou1_B0(0,1,0)*AVG_1_1_LEV_42 + u1_B0(0,1,0)*rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

   CS_03 = rho_B0(0,1,0)*AVG_1_1_LEV_00 + rhou2_B0(0,1,0)*AVG_1_1_LEV_03;

    CS_13 = rho_B0(0,1,0)*AVG_1_1_LEV_10 + rhoE_B0(0,1,0)*AVG_1_1_LEV_14 + rhou0_B0(0,1,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,1,0)*AVG_1_1_LEV_12 + rhou2_B0(0,1,0)*AVG_1_1_LEV_13;

   CS_23 = rho_B0(0,1,0)*AVG_1_1_LEV_20 + rhou0_B0(0,1,0)*AVG_1_1_LEV_21;

    CS_33 = rho_B0(0,1,0)*AVG_1_1_LEV_30 + rhoE_B0(0,1,0)*AVG_1_1_LEV_34 + rhou0_B0(0,1,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,1,0)*AVG_1_1_LEV_32 + rhou2_B0(0,1,0)*AVG_1_1_LEV_33;

    CS_43 = rho_B0(0,1,0)*AVG_1_1_LEV_40 + rhoE_B0(0,1,0)*AVG_1_1_LEV_44 + rhou0_B0(0,1,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,1,0)*AVG_1_1_LEV_42 + rhou2_B0(0,1,0)*AVG_1_1_LEV_43;

   CF_04 = rhou1_B0(0,2,0)*AVG_1_1_LEV_00 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_03;

    CF_14 = p_B0(0,2,0)*AVG_1_1_LEV_12 + rhou1_B0(0,2,0)*AVG_1_1_LEV_10 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_14 +
      u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_14 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_11 +
      u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_12 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

   CF_24 = rhou1_B0(0,2,0)*AVG_1_1_LEV_20 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_21;

    CF_34 = p_B0(0,2,0)*AVG_1_1_LEV_32 + rhou1_B0(0,2,0)*AVG_1_1_LEV_30 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_34 +
      u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_34 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_31 +
      u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_32 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

    CF_44 = p_B0(0,2,0)*AVG_1_1_LEV_42 + rhou1_B0(0,2,0)*AVG_1_1_LEV_40 + p_B0(0,2,0)*u1_B0(0,2,0)*AVG_1_1_LEV_44 +
      u1_B0(0,2,0)*rhoE_B0(0,2,0)*AVG_1_1_LEV_44 + u1_B0(0,2,0)*rhou0_B0(0,2,0)*AVG_1_1_LEV_41 +
      u1_B0(0,2,0)*rhou1_B0(0,2,0)*AVG_1_1_LEV_42 + u1_B0(0,2,0)*rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

   CS_04 = rho_B0(0,2,0)*AVG_1_1_LEV_00 + rhou2_B0(0,2,0)*AVG_1_1_LEV_03;

    CS_14 = rho_B0(0,2,0)*AVG_1_1_LEV_10 + rhoE_B0(0,2,0)*AVG_1_1_LEV_14 + rhou0_B0(0,2,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,2,0)*AVG_1_1_LEV_12 + rhou2_B0(0,2,0)*AVG_1_1_LEV_13;

   CS_24 = rho_B0(0,2,0)*AVG_1_1_LEV_20 + rhou0_B0(0,2,0)*AVG_1_1_LEV_21;

    CS_34 = rho_B0(0,2,0)*AVG_1_1_LEV_30 + rhoE_B0(0,2,0)*AVG_1_1_LEV_34 + rhou0_B0(0,2,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,2,0)*AVG_1_1_LEV_32 + rhou2_B0(0,2,0)*AVG_1_1_LEV_33;

    CS_44 = rho_B0(0,2,0)*AVG_1_1_LEV_40 + rhoE_B0(0,2,0)*AVG_1_1_LEV_44 + rhou0_B0(0,2,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,2,0)*AVG_1_1_LEV_42 + rhou2_B0(0,2,0)*AVG_1_1_LEV_43;

   CF_05 = rhou1_B0(0,3,0)*AVG_1_1_LEV_00 + u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_03;

    CF_15 = p_B0(0,3,0)*AVG_1_1_LEV_12 + rhou1_B0(0,3,0)*AVG_1_1_LEV_10 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_14 +
      u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_14 + u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_11 +
      u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_12 + u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_13;

   CF_25 = rhou1_B0(0,3,0)*AVG_1_1_LEV_20 + u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_21;

    CF_35 = p_B0(0,3,0)*AVG_1_1_LEV_32 + rhou1_B0(0,3,0)*AVG_1_1_LEV_30 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_34 +
      u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_34 + u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_31 +
      u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_32 + u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_33;

    CF_45 = p_B0(0,3,0)*AVG_1_1_LEV_42 + rhou1_B0(0,3,0)*AVG_1_1_LEV_40 + p_B0(0,3,0)*u1_B0(0,3,0)*AVG_1_1_LEV_44 +
      u1_B0(0,3,0)*rhoE_B0(0,3,0)*AVG_1_1_LEV_44 + u1_B0(0,3,0)*rhou0_B0(0,3,0)*AVG_1_1_LEV_41 +
      u1_B0(0,3,0)*rhou1_B0(0,3,0)*AVG_1_1_LEV_42 + u1_B0(0,3,0)*rhou2_B0(0,3,0)*AVG_1_1_LEV_43;

   CS_05 = rho_B0(0,3,0)*AVG_1_1_LEV_00 + rhou2_B0(0,3,0)*AVG_1_1_LEV_03;

    CS_15 = rho_B0(0,3,0)*AVG_1_1_LEV_10 + rhoE_B0(0,3,0)*AVG_1_1_LEV_14 + rhou0_B0(0,3,0)*AVG_1_1_LEV_11 +
      rhou1_B0(0,3,0)*AVG_1_1_LEV_12 + rhou2_B0(0,3,0)*AVG_1_1_LEV_13;

   CS_25 = rho_B0(0,3,0)*AVG_1_1_LEV_20 + rhou0_B0(0,3,0)*AVG_1_1_LEV_21;

    CS_35 = rho_B0(0,3,0)*AVG_1_1_LEV_30 + rhoE_B0(0,3,0)*AVG_1_1_LEV_34 + rhou0_B0(0,3,0)*AVG_1_1_LEV_31 +
      rhou1_B0(0,3,0)*AVG_1_1_LEV_32 + rhou2_B0(0,3,0)*AVG_1_1_LEV_33;

    CS_45 = rho_B0(0,3,0)*AVG_1_1_LEV_40 + rhoE_B0(0,3,0)*AVG_1_1_LEV_44 + rhou0_B0(0,3,0)*AVG_1_1_LEV_41 +
      rhou1_B0(0,3,0)*AVG_1_1_LEV_42 + rhou2_B0(0,3,0)*AVG_1_1_LEV_43;

   max_lambda_00 = shock_filter_control*fmax(fabs(u1_B0(0,0,0)), fabs(u1_B0(0,1,0)));

   max_lambda_11 = max_lambda_00;

   max_lambda_22 = max_lambda_00;

   max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,0,0) + u1_B0(0,0,0)), fabs(a_B0(0,1,0) + u1_B0(0,1,0)));

   max_lambda_44 = shock_filter_control*fmax(fabs(-u1_B0(0,1,0) + a_B0(0,1,0)), fabs(-u1_B0(0,0,0) + a_B0(0,0,0)));

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 +
      CF_03))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 + CF_03))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03) -
      (CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 +
      CF_03) - (CS_02*max_lambda_00 + CF_02)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 + CF_03) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 +
      CF_03) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(CS_02*max_lambda_00 +
      CF_02) + ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 + CF_01) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 +
      CF_01) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) +
      ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01))*(((1.0/2.0))*(CS_00*max_lambda_00 +
      CF_00) + ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01)));

    beta_3 = -(781.0/480.0)*(CS_04*max_lambda_00 + CF_04) - (781.0/1440.0)*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/36.0))*((9*(CS_03*max_lambda_00 + CF_03) - (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) -
      (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05)*(9*(CS_03*max_lambda_00 + CF_03) -
      (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) - (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05))
      + ((13.0/12.0))*((2*(CS_04*max_lambda_00 + CF_04) - (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) -
      (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02)*(2*(CS_04*max_lambda_00 + CF_04) -
      (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) - (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02))
      + ((781.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((781.0/1440.0))*(CS_05*max_lambda_00 + CF_05);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(CS_03*max_lambda_00 + CF_03))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/8.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/24.0))*(CS_05*max_lambda_00 + CF_05) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_00*max_lambda_00 + CF_00) +
      ((11.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 + CF_04))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04) -
      (-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
      CF_04) - (-CS_03*max_lambda_00 + CF_03)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00 + CF_02) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00
      + CF_02) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (-CS_02*max_lambda_00 +
      CF_02))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) -
      (-CS_02*max_lambda_00 + CF_02)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00 + CF_04) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00
      + CF_04) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 +
      CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - (-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) -
      (-CS_04*max_lambda_00 + CF_04)));

    beta_3 = ((1.0/36.0))*((-(-CS_00*max_lambda_00 + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) +
      ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(-(-CS_00*max_lambda_00
      + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) + ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) +
      ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((2*(-CS_01*max_lambda_00 + CF_01) -
      (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 +
      CF_03)*(2*(-CS_01*max_lambda_00 + CF_01) - (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 + CF_03)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 + CF_01)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(-CS_02*max_lambda_00 + CF_02))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((1.0/24.0))*(-CS_00*max_lambda_00 + CF_00) + ((13.0/24.0))*(-CS_02*max_lambda_00 + CF_02))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_05*max_lambda_00 + CF_05) +
      ((11.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 +
      CF_13))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 + CF_13))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13) -
      (CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 +
      CF_13) - (CS_12*max_lambda_11 + CF_12)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 + CF_13) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 +
      CF_13) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(CS_12*max_lambda_11 +
      CF_12) + ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 + CF_11) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 +
      CF_11) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) +
      ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11))*(((1.0/2.0))*(CS_10*max_lambda_11 +
      CF_10) + ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11)));

    beta_3 = -(781.0/480.0)*(CS_14*max_lambda_11 + CF_14) - (781.0/1440.0)*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/36.0))*((9*(CS_13*max_lambda_11 + CF_13) - (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) -
      (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15)*(9*(CS_13*max_lambda_11 + CF_13) -
      (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) - (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15))
      + ((13.0/12.0))*((2*(CS_14*max_lambda_11 + CF_14) - (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) -
      (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12)*(2*(CS_14*max_lambda_11 + CF_14) -
      (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) - (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12))
      + ((781.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((781.0/1440.0))*(CS_15*max_lambda_11 + CF_15);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(CS_13*max_lambda_11 + CF_13))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/8.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/24.0))*(CS_15*max_lambda_11 + CF_15) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_10*max_lambda_11 + CF_10) +
      ((11.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 + CF_14))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14) -
      (-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
      CF_14) - (-CS_13*max_lambda_11 + CF_13)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11 + CF_12) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11
      + CF_12) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (-CS_12*max_lambda_11 +
      CF_12))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) -
      (-CS_12*max_lambda_11 + CF_12)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11 + CF_14) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11
      + CF_14) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 +
      CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - (-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) -
      (-CS_14*max_lambda_11 + CF_14)));

    beta_3 = ((1.0/36.0))*((-(-CS_10*max_lambda_11 + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) +
      ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(-(-CS_10*max_lambda_11
      + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) + ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) +
      ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((2*(-CS_11*max_lambda_11 + CF_11) -
      (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 +
      CF_13)*(2*(-CS_11*max_lambda_11 + CF_11) - (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 + CF_13)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 + CF_11)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(-CS_12*max_lambda_11 + CF_12))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((1.0/24.0))*(-CS_10*max_lambda_11 + CF_10) + ((13.0/24.0))*(-CS_12*max_lambda_11 + CF_12))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_15*max_lambda_11 + CF_15) +
      ((11.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 +
      CF_23))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 + CF_23))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23) -
      (CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 +
      CF_23) - (CS_22*max_lambda_22 + CF_22)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 + CF_23) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 +
      CF_23) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(CS_22*max_lambda_22 +
      CF_22) + ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 + CF_21) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 +
      CF_21) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) +
      ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21))*(((1.0/2.0))*(CS_20*max_lambda_22 +
      CF_20) + ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21)));

    beta_3 = -(781.0/480.0)*(CS_24*max_lambda_22 + CF_24) - (781.0/1440.0)*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/36.0))*((9*(CS_23*max_lambda_22 + CF_23) - (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) -
      (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25)*(9*(CS_23*max_lambda_22 + CF_23) -
      (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) - (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25))
      + ((13.0/12.0))*((2*(CS_24*max_lambda_22 + CF_24) - (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) -
      (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22)*(2*(CS_24*max_lambda_22 + CF_24) -
      (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) - (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22))
      + ((781.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((781.0/1440.0))*(CS_25*max_lambda_22 + CF_25);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(CS_23*max_lambda_22 + CF_23))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/8.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/24.0))*(CS_25*max_lambda_22 + CF_25) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_20*max_lambda_22 + CF_20) +
      ((11.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 + CF_24))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24) -
      (-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
      CF_24) - (-CS_23*max_lambda_22 + CF_23)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22 + CF_22) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22
      + CF_22) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (-CS_22*max_lambda_22 +
      CF_22))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) -
      (-CS_22*max_lambda_22 + CF_22)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22 + CF_24) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22
      + CF_24) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 +
      CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - (-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) -
      (-CS_24*max_lambda_22 + CF_24)));

    beta_3 = ((1.0/36.0))*((-(-CS_20*max_lambda_22 + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) +
      ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(-(-CS_20*max_lambda_22
      + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) + ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) +
      ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((2*(-CS_21*max_lambda_22 + CF_21) -
      (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 +
      CF_23)*(2*(-CS_21*max_lambda_22 + CF_21) - (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 + CF_23)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 + CF_21)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(-CS_22*max_lambda_22 + CF_22))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((1.0/24.0))*(-CS_20*max_lambda_22 + CF_20) + ((13.0/24.0))*(-CS_22*max_lambda_22 + CF_22))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_25*max_lambda_22 + CF_25) +
      ((11.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 +
      CF_33))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 + CF_33))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33) -
      (CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 +
      CF_33) - (CS_32*max_lambda_33 + CF_32)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 + CF_33) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 +
      CF_33) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(CS_32*max_lambda_33 +
      CF_32) + ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 + CF_31) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 +
      CF_31) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) +
      ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31))*(((1.0/2.0))*(CS_30*max_lambda_33 +
      CF_30) + ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31)));

    beta_3 = -(781.0/480.0)*(CS_34*max_lambda_33 + CF_34) - (781.0/1440.0)*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/36.0))*((9*(CS_33*max_lambda_33 + CF_33) - (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) -
      (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35)*(9*(CS_33*max_lambda_33 + CF_33) -
      (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) - (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35))
      + ((13.0/12.0))*((2*(CS_34*max_lambda_33 + CF_34) - (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) -
      (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32)*(2*(CS_34*max_lambda_33 + CF_34) -
      (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) - (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32))
      + ((781.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((781.0/1440.0))*(CS_35*max_lambda_33 + CF_35);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(CS_33*max_lambda_33 + CF_33))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/8.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/24.0))*(CS_35*max_lambda_33 + CF_35) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_30*max_lambda_33 + CF_30) +
      ((11.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 + CF_34))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34) -
      (-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
      CF_34) - (-CS_33*max_lambda_33 + CF_33)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33 + CF_32) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33
      + CF_32) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (-CS_32*max_lambda_33 +
      CF_32))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) -
      (-CS_32*max_lambda_33 + CF_32)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33 + CF_34) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33
      + CF_34) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 +
      CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - (-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) -
      (-CS_34*max_lambda_33 + CF_34)));

    beta_3 = ((1.0/36.0))*((-(-CS_30*max_lambda_33 + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) +
      ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(-(-CS_30*max_lambda_33
      + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) + ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) +
      ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((2*(-CS_31*max_lambda_33 + CF_31) -
      (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 +
      CF_33)*(2*(-CS_31*max_lambda_33 + CF_31) - (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 + CF_33)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 + CF_31)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(-CS_32*max_lambda_33 + CF_32))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((1.0/24.0))*(-CS_30*max_lambda_33 + CF_30) + ((13.0/24.0))*(-CS_32*max_lambda_33 + CF_32))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_35*max_lambda_33 + CF_35) +
      ((11.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 +
      CF_43))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 + CF_43))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43) -
      (CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 +
      CF_43) - (CS_42*max_lambda_44 + CF_42)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 + CF_43) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 +
      CF_43) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(CS_42*max_lambda_44 +
      CF_42) + ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 + CF_41) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 +
      CF_41) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) +
      ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41))*(((1.0/2.0))*(CS_40*max_lambda_44 +
      CF_40) + ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41)));

    beta_3 = -(781.0/480.0)*(CS_44*max_lambda_44 + CF_44) - (781.0/1440.0)*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/36.0))*((9*(CS_43*max_lambda_44 + CF_43) - (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) -
      (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45)*(9*(CS_43*max_lambda_44 + CF_43) -
      (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) - (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45))
      + ((13.0/12.0))*((2*(CS_44*max_lambda_44 + CF_44) - (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) -
      (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42)*(2*(CS_44*max_lambda_44 + CF_44) -
      (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) - (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42))
      + ((781.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((781.0/1440.0))*(CS_45*max_lambda_44 + CF_45);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(CS_43*max_lambda_44 + CF_43))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/8.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/24.0))*(CS_45*max_lambda_44 + CF_45) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_40*max_lambda_44 + CF_40) +
      ((11.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_0*inv_omega_sum + Recon_4;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 + CF_44))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44) -
      (-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
      CF_44) - (-CS_43*max_lambda_44 + CF_43)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44 + CF_42) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44
      + CF_42) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (-CS_42*max_lambda_44 +
      CF_42))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) -
      (-CS_42*max_lambda_44 + CF_42)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44 + CF_44) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44
      + CF_44) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 +
      CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - (-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) -
      (-CS_44*max_lambda_44 + CF_44)));

    beta_3 = ((1.0/36.0))*((-(-CS_40*max_lambda_44 + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) +
      ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(-(-CS_40*max_lambda_44
      + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) + ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) +
      ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((2*(-CS_41*max_lambda_44 + CF_41) -
      (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 +
      CF_43)*(2*(-CS_41*max_lambda_44 + CF_41) - (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 + CF_43)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 + CF_41)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(-CS_42*max_lambda_44 + CF_42))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((1.0/24.0))*(-CS_40*max_lambda_44 + CF_40) + ((13.0/24.0))*(-CS_42*max_lambda_44 + CF_42))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_45*max_lambda_44 + CF_45) +
      ((11.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_0*inv_omega_sum + Recon_4;

    wk5_B0(0,0,0) = 0.707106781186547*AVG_1_rho*Recon_3*inv_AVG_a + 0.707106781186547*AVG_1_rho*Recon_4*inv_AVG_a +
      Recon_1;

    wk6_B0(0,0,0) = AVG_1_rho*Recon_2 + AVG_1_u0*Recon_1 + 0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_3*inv_AVG_a +
      0.707106781186547*AVG_1_rho*AVG_1_u0*Recon_4*inv_AVG_a;

    wk7_B0(0,0,0) = AVG_1_u1*Recon_1 + 0.707106781186547*(-AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a +
      0.707106781186547*(AVG_1_a + AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a;

    wk8_B0(0,0,0) = AVG_1_u2*Recon_1 - AVG_1_rho*Recon_0 + 0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_3*inv_AVG_a +
      0.707106781186547*AVG_1_rho*AVG_1_u2*Recon_4*inv_AVG_a;

    wk9_B0(0,0,0) = (((1.0/2.0))*(AVG_1_u0*AVG_1_u0) + ((1.0/2.0))*(AVG_1_u1*AVG_1_u1) +
      ((1.0/2.0))*(AVG_1_u2*AVG_1_u2))*Recon_1 + AVG_1_rho*AVG_1_u0*Recon_2 - AVG_1_rho*AVG_1_u2*Recon_0 +
      0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
      (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 + AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_3*inv_AVG_a +
      0.707106781186547*(((AVG_1_a*AVG_1_a) + ((1.0/2.0))*((AVG_1_u0*AVG_1_u0) + (AVG_1_u1*AVG_1_u1) +
      (AVG_1_u2*AVG_1_u2))*gamma_m1)*invgamma_m1 - AVG_1_a*AVG_1_u1)*AVG_1_rho*Recon_4*inv_AVG_a;

}

 void opensbliblock00Kernel002(const ACC<double> &a_B0, const ACC<double> &p_B0, const ACC<double> &rhoE_B0, const
ACC<double> &rho_B0, const ACC<double> &rhou0_B0, const ACC<double> &rhou1_B0, const ACC<double> &rhou2_B0, const
ACC<double> &u0_B0, const ACC<double> &u1_B0, const ACC<double> &u2_B0, ACC<double> &wk10_B0, ACC<double> &wk11_B0,
ACC<double> &wk12_B0, ACC<double> &wk13_B0, ACC<double> &wk14_B0)
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
   double CF_04 = 0.0;
   double CF_05 = 0.0;
   double CF_10 = 0.0;
   double CF_11 = 0.0;
   double CF_12 = 0.0;
   double CF_13 = 0.0;
   double CF_14 = 0.0;
   double CF_15 = 0.0;
   double CF_20 = 0.0;
   double CF_21 = 0.0;
   double CF_22 = 0.0;
   double CF_23 = 0.0;
   double CF_24 = 0.0;
   double CF_25 = 0.0;
   double CF_30 = 0.0;
   double CF_31 = 0.0;
   double CF_32 = 0.0;
   double CF_33 = 0.0;
   double CF_34 = 0.0;
   double CF_35 = 0.0;
   double CF_40 = 0.0;
   double CF_41 = 0.0;
   double CF_42 = 0.0;
   double CF_43 = 0.0;
   double CF_44 = 0.0;
   double CF_45 = 0.0;
   double CS_00 = 0.0;
   double CS_01 = 0.0;
   double CS_02 = 0.0;
   double CS_03 = 0.0;
   double CS_04 = 0.0;
   double CS_05 = 0.0;
   double CS_10 = 0.0;
   double CS_11 = 0.0;
   double CS_12 = 0.0;
   double CS_13 = 0.0;
   double CS_14 = 0.0;
   double CS_15 = 0.0;
   double CS_20 = 0.0;
   double CS_21 = 0.0;
   double CS_22 = 0.0;
   double CS_23 = 0.0;
   double CS_24 = 0.0;
   double CS_25 = 0.0;
   double CS_30 = 0.0;
   double CS_31 = 0.0;
   double CS_32 = 0.0;
   double CS_33 = 0.0;
   double CS_34 = 0.0;
   double CS_35 = 0.0;
   double CS_40 = 0.0;
   double CS_41 = 0.0;
   double CS_42 = 0.0;
   double CS_43 = 0.0;
   double CS_44 = 0.0;
   double CS_45 = 0.0;
   double Recon_0 = 0.0;
   double Recon_1 = 0.0;
   double Recon_2 = 0.0;
   double Recon_3 = 0.0;
   double Recon_4 = 0.0;
   double alpha_0 = 0.0;
   double alpha_1 = 0.0;
   double alpha_2 = 0.0;
   double alpha_3 = 0.0;
   double beta_0 = 0.0;
   double beta_1 = 0.0;
   double beta_2 = 0.0;
   double beta_3 = 0.0;
   double delta_0 = 0.0;
   double delta_1 = 0.0;
   double delta_2 = 0.0;
   double delta_3 = 0.0;
   double inv_AVG_a = 0.0;
   double inv_AVG_rho = 0.0;
   double inv_alpha_sum = 0.0;
   double inv_beta_0 = 0.0;
   double inv_beta_1 = 0.0;
   double inv_beta_2 = 0.0;
   double inv_beta_3 = 0.0;
   double inv_omega_sum = 0.0;
   double max_lambda_00 = 0.0;
   double max_lambda_11 = 0.0;
   double max_lambda_22 = 0.0;
   double max_lambda_33 = 0.0;
   double max_lambda_44 = 0.0;
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

   CF_00 = rhou2_B0(0,0,-2)*AVG_2_2_LEV_00 + u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_02;

   CF_10 = rhou2_B0(0,0,-2)*AVG_2_2_LEV_10 + u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_11;

    CF_20 = p_B0(0,0,-2)*AVG_2_2_LEV_23 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_20 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_24 +
      u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_24 + u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_21 +
      u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_22 + u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_23;

    CF_30 = p_B0(0,0,-2)*AVG_2_2_LEV_33 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_30 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_34 +
      u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_34 + u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_31 +
      u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_32 + u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_33;

    CF_40 = p_B0(0,0,-2)*AVG_2_2_LEV_43 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_40 + p_B0(0,0,-2)*u2_B0(0,0,-2)*AVG_2_2_LEV_44 +
      u2_B0(0,0,-2)*rhoE_B0(0,0,-2)*AVG_2_2_LEV_44 + u2_B0(0,0,-2)*rhou0_B0(0,0,-2)*AVG_2_2_LEV_41 +
      u2_B0(0,0,-2)*rhou1_B0(0,0,-2)*AVG_2_2_LEV_42 + u2_B0(0,0,-2)*rhou2_B0(0,0,-2)*AVG_2_2_LEV_43;

   CS_00 = rho_B0(0,0,-2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-2)*AVG_2_2_LEV_02;

   CS_10 = rho_B0(0,0,-2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_11;

    CS_20 = rho_B0(0,0,-2)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,-2)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_23;

    CS_30 = rho_B0(0,0,-2)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,-2)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_33;

    CS_40 = rho_B0(0,0,-2)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-2)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-2)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,-2)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-2)*AVG_2_2_LEV_43;

   CF_01 = rhou2_B0(0,0,-1)*AVG_2_2_LEV_00 + u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_02;

   CF_11 = rhou2_B0(0,0,-1)*AVG_2_2_LEV_10 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_11;

    CF_21 = p_B0(0,0,-1)*AVG_2_2_LEV_23 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_20 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_24 +
      u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 +
      u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

    CF_31 = p_B0(0,0,-1)*AVG_2_2_LEV_33 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_30 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_34 +
      u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 +
      u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

    CF_41 = p_B0(0,0,-1)*AVG_2_2_LEV_43 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_40 + p_B0(0,0,-1)*u2_B0(0,0,-1)*AVG_2_2_LEV_44 +
      u2_B0(0,0,-1)*rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 + u2_B0(0,0,-1)*rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 +
      u2_B0(0,0,-1)*rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 + u2_B0(0,0,-1)*rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

   CS_01 = rho_B0(0,0,-1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,-1)*AVG_2_2_LEV_02;

   CS_11 = rho_B0(0,0,-1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_11;

    CS_21 = rho_B0(0,0,-1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,-1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_23;

    CS_31 = rho_B0(0,0,-1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,-1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_33;

    CS_41 = rho_B0(0,0,-1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,-1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,-1)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,-1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,-1)*AVG_2_2_LEV_43;

   CF_02 = rhou2_B0(0,0,0)*AVG_2_2_LEV_00 + u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_02;

   CF_12 = rhou2_B0(0,0,0)*AVG_2_2_LEV_10 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_11;

    CF_22 = p_B0(0,0,0)*AVG_2_2_LEV_23 + rhou2_B0(0,0,0)*AVG_2_2_LEV_20 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_24 +
      u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_24 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_21 +
      u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_22 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

    CF_32 = p_B0(0,0,0)*AVG_2_2_LEV_33 + rhou2_B0(0,0,0)*AVG_2_2_LEV_30 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_34 +
      u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_34 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_31 +
      u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_32 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

    CF_42 = p_B0(0,0,0)*AVG_2_2_LEV_43 + rhou2_B0(0,0,0)*AVG_2_2_LEV_40 + p_B0(0,0,0)*u2_B0(0,0,0)*AVG_2_2_LEV_44 +
      u2_B0(0,0,0)*rhoE_B0(0,0,0)*AVG_2_2_LEV_44 + u2_B0(0,0,0)*rhou0_B0(0,0,0)*AVG_2_2_LEV_41 +
      u2_B0(0,0,0)*rhou1_B0(0,0,0)*AVG_2_2_LEV_42 + u2_B0(0,0,0)*rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

   CS_02 = rho_B0(0,0,0)*AVG_2_2_LEV_00 + rhou1_B0(0,0,0)*AVG_2_2_LEV_02;

   CS_12 = rho_B0(0,0,0)*AVG_2_2_LEV_10 + rhou0_B0(0,0,0)*AVG_2_2_LEV_11;

    CS_22 = rho_B0(0,0,0)*AVG_2_2_LEV_20 + rhoE_B0(0,0,0)*AVG_2_2_LEV_24 + rhou0_B0(0,0,0)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,0)*AVG_2_2_LEV_22 + rhou2_B0(0,0,0)*AVG_2_2_LEV_23;

    CS_32 = rho_B0(0,0,0)*AVG_2_2_LEV_30 + rhoE_B0(0,0,0)*AVG_2_2_LEV_34 + rhou0_B0(0,0,0)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,0)*AVG_2_2_LEV_32 + rhou2_B0(0,0,0)*AVG_2_2_LEV_33;

    CS_42 = rho_B0(0,0,0)*AVG_2_2_LEV_40 + rhoE_B0(0,0,0)*AVG_2_2_LEV_44 + rhou0_B0(0,0,0)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,0)*AVG_2_2_LEV_42 + rhou2_B0(0,0,0)*AVG_2_2_LEV_43;

   CF_03 = rhou2_B0(0,0,1)*AVG_2_2_LEV_00 + u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_02;

   CF_13 = rhou2_B0(0,0,1)*AVG_2_2_LEV_10 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_11;

    CF_23 = p_B0(0,0,1)*AVG_2_2_LEV_23 + rhou2_B0(0,0,1)*AVG_2_2_LEV_20 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_24 +
      u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_24 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_21 +
      u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_22 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

    CF_33 = p_B0(0,0,1)*AVG_2_2_LEV_33 + rhou2_B0(0,0,1)*AVG_2_2_LEV_30 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_34 +
      u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_34 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_31 +
      u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_32 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

    CF_43 = p_B0(0,0,1)*AVG_2_2_LEV_43 + rhou2_B0(0,0,1)*AVG_2_2_LEV_40 + p_B0(0,0,1)*u2_B0(0,0,1)*AVG_2_2_LEV_44 +
      u2_B0(0,0,1)*rhoE_B0(0,0,1)*AVG_2_2_LEV_44 + u2_B0(0,0,1)*rhou0_B0(0,0,1)*AVG_2_2_LEV_41 +
      u2_B0(0,0,1)*rhou1_B0(0,0,1)*AVG_2_2_LEV_42 + u2_B0(0,0,1)*rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

   CS_03 = rho_B0(0,0,1)*AVG_2_2_LEV_00 + rhou1_B0(0,0,1)*AVG_2_2_LEV_02;

   CS_13 = rho_B0(0,0,1)*AVG_2_2_LEV_10 + rhou0_B0(0,0,1)*AVG_2_2_LEV_11;

    CS_23 = rho_B0(0,0,1)*AVG_2_2_LEV_20 + rhoE_B0(0,0,1)*AVG_2_2_LEV_24 + rhou0_B0(0,0,1)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,1)*AVG_2_2_LEV_22 + rhou2_B0(0,0,1)*AVG_2_2_LEV_23;

    CS_33 = rho_B0(0,0,1)*AVG_2_2_LEV_30 + rhoE_B0(0,0,1)*AVG_2_2_LEV_34 + rhou0_B0(0,0,1)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,1)*AVG_2_2_LEV_32 + rhou2_B0(0,0,1)*AVG_2_2_LEV_33;

    CS_43 = rho_B0(0,0,1)*AVG_2_2_LEV_40 + rhoE_B0(0,0,1)*AVG_2_2_LEV_44 + rhou0_B0(0,0,1)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,1)*AVG_2_2_LEV_42 + rhou2_B0(0,0,1)*AVG_2_2_LEV_43;

   CF_04 = rhou2_B0(0,0,2)*AVG_2_2_LEV_00 + u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_02;

   CF_14 = rhou2_B0(0,0,2)*AVG_2_2_LEV_10 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_11;

    CF_24 = p_B0(0,0,2)*AVG_2_2_LEV_23 + rhou2_B0(0,0,2)*AVG_2_2_LEV_20 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_24 +
      u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_24 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_21 +
      u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_22 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

    CF_34 = p_B0(0,0,2)*AVG_2_2_LEV_33 + rhou2_B0(0,0,2)*AVG_2_2_LEV_30 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_34 +
      u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_34 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_31 +
      u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_32 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

    CF_44 = p_B0(0,0,2)*AVG_2_2_LEV_43 + rhou2_B0(0,0,2)*AVG_2_2_LEV_40 + p_B0(0,0,2)*u2_B0(0,0,2)*AVG_2_2_LEV_44 +
      u2_B0(0,0,2)*rhoE_B0(0,0,2)*AVG_2_2_LEV_44 + u2_B0(0,0,2)*rhou0_B0(0,0,2)*AVG_2_2_LEV_41 +
      u2_B0(0,0,2)*rhou1_B0(0,0,2)*AVG_2_2_LEV_42 + u2_B0(0,0,2)*rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

   CS_04 = rho_B0(0,0,2)*AVG_2_2_LEV_00 + rhou1_B0(0,0,2)*AVG_2_2_LEV_02;

   CS_14 = rho_B0(0,0,2)*AVG_2_2_LEV_10 + rhou0_B0(0,0,2)*AVG_2_2_LEV_11;

    CS_24 = rho_B0(0,0,2)*AVG_2_2_LEV_20 + rhoE_B0(0,0,2)*AVG_2_2_LEV_24 + rhou0_B0(0,0,2)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,2)*AVG_2_2_LEV_22 + rhou2_B0(0,0,2)*AVG_2_2_LEV_23;

    CS_34 = rho_B0(0,0,2)*AVG_2_2_LEV_30 + rhoE_B0(0,0,2)*AVG_2_2_LEV_34 + rhou0_B0(0,0,2)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,2)*AVG_2_2_LEV_32 + rhou2_B0(0,0,2)*AVG_2_2_LEV_33;

    CS_44 = rho_B0(0,0,2)*AVG_2_2_LEV_40 + rhoE_B0(0,0,2)*AVG_2_2_LEV_44 + rhou0_B0(0,0,2)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,2)*AVG_2_2_LEV_42 + rhou2_B0(0,0,2)*AVG_2_2_LEV_43;

   CF_05 = rhou2_B0(0,0,3)*AVG_2_2_LEV_00 + u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_02;

   CF_15 = rhou2_B0(0,0,3)*AVG_2_2_LEV_10 + u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_11;

    CF_25 = p_B0(0,0,3)*AVG_2_2_LEV_23 + rhou2_B0(0,0,3)*AVG_2_2_LEV_20 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_24 +
      u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_24 + u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_21 +
      u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_22 + u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_23;

    CF_35 = p_B0(0,0,3)*AVG_2_2_LEV_33 + rhou2_B0(0,0,3)*AVG_2_2_LEV_30 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_34 +
      u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_34 + u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_31 +
      u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_32 + u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_33;

    CF_45 = p_B0(0,0,3)*AVG_2_2_LEV_43 + rhou2_B0(0,0,3)*AVG_2_2_LEV_40 + p_B0(0,0,3)*u2_B0(0,0,3)*AVG_2_2_LEV_44 +
      u2_B0(0,0,3)*rhoE_B0(0,0,3)*AVG_2_2_LEV_44 + u2_B0(0,0,3)*rhou0_B0(0,0,3)*AVG_2_2_LEV_41 +
      u2_B0(0,0,3)*rhou1_B0(0,0,3)*AVG_2_2_LEV_42 + u2_B0(0,0,3)*rhou2_B0(0,0,3)*AVG_2_2_LEV_43;

   CS_05 = rho_B0(0,0,3)*AVG_2_2_LEV_00 + rhou1_B0(0,0,3)*AVG_2_2_LEV_02;

   CS_15 = rho_B0(0,0,3)*AVG_2_2_LEV_10 + rhou0_B0(0,0,3)*AVG_2_2_LEV_11;

    CS_25 = rho_B0(0,0,3)*AVG_2_2_LEV_20 + rhoE_B0(0,0,3)*AVG_2_2_LEV_24 + rhou0_B0(0,0,3)*AVG_2_2_LEV_21 +
      rhou1_B0(0,0,3)*AVG_2_2_LEV_22 + rhou2_B0(0,0,3)*AVG_2_2_LEV_23;

    CS_35 = rho_B0(0,0,3)*AVG_2_2_LEV_30 + rhoE_B0(0,0,3)*AVG_2_2_LEV_34 + rhou0_B0(0,0,3)*AVG_2_2_LEV_31 +
      rhou1_B0(0,0,3)*AVG_2_2_LEV_32 + rhou2_B0(0,0,3)*AVG_2_2_LEV_33;

    CS_45 = rho_B0(0,0,3)*AVG_2_2_LEV_40 + rhoE_B0(0,0,3)*AVG_2_2_LEV_44 + rhou0_B0(0,0,3)*AVG_2_2_LEV_41 +
      rhou1_B0(0,0,3)*AVG_2_2_LEV_42 + rhou2_B0(0,0,3)*AVG_2_2_LEV_43;

   max_lambda_00 = shock_filter_control*fmax(fabs(u2_B0(0,0,1)), fabs(u2_B0(0,0,0)));

   max_lambda_11 = max_lambda_00;

   max_lambda_22 = max_lambda_00;

   max_lambda_33 = shock_filter_control*fmax(fabs(a_B0(0,0,1) + u2_B0(0,0,1)), fabs(a_B0(0,0,0) + u2_B0(0,0,0)));

   max_lambda_44 = shock_filter_control*fmax(fabs(-u2_B0(0,0,0) + a_B0(0,0,0)), fabs(-u2_B0(0,0,1) + a_B0(0,0,1)));

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 +
      CF_03))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) - (1.0/2.0)*(CS_03*max_lambda_00 + CF_03))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 + CF_03) -
      (CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(CS_03*max_lambda_00 +
      CF_03) - (CS_02*max_lambda_00 + CF_02)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 + CF_03) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - 2*(CS_03*max_lambda_00 +
      CF_03) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(CS_02*max_lambda_00 +
      CF_02) + ((1.0/2.0))*(CS_04*max_lambda_00 + CF_04) - (CS_03*max_lambda_00 + CF_03)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 + CF_01) +
      ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))*(((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) - 2*(CS_01*max_lambda_00 +
      CF_01) + ((3.0/2.0))*(CS_02*max_lambda_00 + CF_02))) + ((13.0/12.0))*((((1.0/2.0))*(CS_00*max_lambda_00 + CF_00) +
      ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01))*(((1.0/2.0))*(CS_00*max_lambda_00 +
      CF_00) + ((1.0/2.0))*(CS_02*max_lambda_00 + CF_02) - (CS_01*max_lambda_00 + CF_01)));

    beta_3 = -(781.0/480.0)*(CS_04*max_lambda_00 + CF_04) - (781.0/1440.0)*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/36.0))*((9*(CS_03*max_lambda_00 + CF_03) - (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) -
      (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05)*(9*(CS_03*max_lambda_00 + CF_03) -
      (11.0/2.0)*(CS_02*max_lambda_00 + CF_02) - (9.0/2.0)*(CS_04*max_lambda_00 + CF_04) + CS_05*max_lambda_00 + CF_05))
      + ((13.0/12.0))*((2*(CS_04*max_lambda_00 + CF_04) - (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) -
      (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02)*(2*(CS_04*max_lambda_00 + CF_04) -
      (5.0/2.0)*(CS_03*max_lambda_00 + CF_03) - (1.0/2.0)*(CS_05*max_lambda_00 + CF_05) + CS_02*max_lambda_00 + CF_02))
      + ((781.0/480.0))*(CS_03*max_lambda_00 + CF_03) + ((781.0/1440.0))*(CS_05*max_lambda_00 + CF_05);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(CS_03*max_lambda_00 + CF_03))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_04*max_lambda_00 + CF_04) + ((1.0/8.0))*(CS_02*max_lambda_00 + CF_02) +
      ((1.0/24.0))*(CS_05*max_lambda_00 + CF_05) + ((13.0/24.0))*(CS_03*max_lambda_00 + CF_03))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_00*max_lambda_00 + CF_00) +
      ((11.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(CS_02*max_lambda_00 + CF_02))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_04*max_lambda_00 + CF_04))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 + CF_04) -
      (-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_02*max_lambda_00 + CF_02) + ((1.0/2.0))*(-CS_04*max_lambda_00 +
      CF_04) - (-CS_03*max_lambda_00 + CF_03)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00 + CF_02) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) - 2*(-CS_02*max_lambda_00
      + CF_02) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (-CS_02*max_lambda_00 +
      CF_02))*(((1.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) -
      (-CS_02*max_lambda_00 + CF_02)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00 + CF_04) +
      ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - 2*(-CS_04*max_lambda_00
      + CF_04) + ((3.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 +
      CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) - (-CS_04*max_lambda_00 +
      CF_04))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) + ((1.0/2.0))*(-CS_05*max_lambda_00 + CF_05) -
      (-CS_04*max_lambda_00 + CF_04)));

    beta_3 = ((1.0/36.0))*((-(-CS_00*max_lambda_00 + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) +
      ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) + ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))*(-(-CS_00*max_lambda_00
      + CF_00) - 9*(-CS_02*max_lambda_00 + CF_02) + ((9.0/2.0))*(-CS_01*max_lambda_00 + CF_01) +
      ((11.0/2.0))*(-CS_03*max_lambda_00 + CF_03))) + ((13.0/12.0))*((2*(-CS_01*max_lambda_00 + CF_01) -
      (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) - (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 +
      CF_03)*(2*(-CS_01*max_lambda_00 + CF_01) - (5.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) - CS_03*max_lambda_00 + CF_03)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 +
      CF_01))*(((1.0/2.0))*(-CS_03*max_lambda_00 + CF_03) - (3.0/2.0)*(-CS_02*max_lambda_00 + CF_02) -
      (1.0/2.0)*(-CS_00*max_lambda_00 + CF_00) + ((3.0/2.0))*(-CS_01*max_lambda_00 + CF_01)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_0 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/6.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((5.0/12.0))*(-CS_02*max_lambda_00 + CF_02))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_01*max_lambda_00 + CF_01) + ((1.0/8.0))*(-CS_03*max_lambda_00 + CF_03) +
      ((1.0/24.0))*(-CS_00*max_lambda_00 + CF_00) + ((13.0/24.0))*(-CS_02*max_lambda_00 + CF_02))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_05*max_lambda_00 + CF_05) +
      ((11.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_04*max_lambda_00 + CF_04) + ((1.0/6.0))*(-CS_02*max_lambda_00 + CF_02) +
      ((5.0/12.0))*(-CS_03*max_lambda_00 + CF_03))*delta_0*inv_omega_sum + Recon_0;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 +
      CF_13))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) - (1.0/2.0)*(CS_13*max_lambda_11 + CF_13))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 + CF_13) -
      (CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(CS_13*max_lambda_11 +
      CF_13) - (CS_12*max_lambda_11 + CF_12)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 + CF_13) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - 2*(CS_13*max_lambda_11 +
      CF_13) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(CS_12*max_lambda_11 +
      CF_12) + ((1.0/2.0))*(CS_14*max_lambda_11 + CF_14) - (CS_13*max_lambda_11 + CF_13)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 + CF_11) +
      ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))*(((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) - 2*(CS_11*max_lambda_11 +
      CF_11) + ((3.0/2.0))*(CS_12*max_lambda_11 + CF_12))) + ((13.0/12.0))*((((1.0/2.0))*(CS_10*max_lambda_11 + CF_10) +
      ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11))*(((1.0/2.0))*(CS_10*max_lambda_11 +
      CF_10) + ((1.0/2.0))*(CS_12*max_lambda_11 + CF_12) - (CS_11*max_lambda_11 + CF_11)));

    beta_3 = -(781.0/480.0)*(CS_14*max_lambda_11 + CF_14) - (781.0/1440.0)*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/36.0))*((9*(CS_13*max_lambda_11 + CF_13) - (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) -
      (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15)*(9*(CS_13*max_lambda_11 + CF_13) -
      (11.0/2.0)*(CS_12*max_lambda_11 + CF_12) - (9.0/2.0)*(CS_14*max_lambda_11 + CF_14) + CS_15*max_lambda_11 + CF_15))
      + ((13.0/12.0))*((2*(CS_14*max_lambda_11 + CF_14) - (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) -
      (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12)*(2*(CS_14*max_lambda_11 + CF_14) -
      (5.0/2.0)*(CS_13*max_lambda_11 + CF_13) - (1.0/2.0)*(CS_15*max_lambda_11 + CF_15) + CS_12*max_lambda_11 + CF_12))
      + ((781.0/480.0))*(CS_13*max_lambda_11 + CF_13) + ((781.0/1440.0))*(CS_15*max_lambda_11 + CF_15);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(CS_13*max_lambda_11 + CF_13))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_14*max_lambda_11 + CF_14) + ((1.0/8.0))*(CS_12*max_lambda_11 + CF_12) +
      ((1.0/24.0))*(CS_15*max_lambda_11 + CF_15) + ((13.0/24.0))*(CS_13*max_lambda_11 + CF_13))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_10*max_lambda_11 + CF_10) +
      ((11.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(CS_12*max_lambda_11 + CF_12))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_14*max_lambda_11 + CF_14))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 + CF_14) -
      (-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_12*max_lambda_11 + CF_12) + ((1.0/2.0))*(-CS_14*max_lambda_11 +
      CF_14) - (-CS_13*max_lambda_11 + CF_13)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11 + CF_12) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) - 2*(-CS_12*max_lambda_11
      + CF_12) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (-CS_12*max_lambda_11 +
      CF_12))*(((1.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) -
      (-CS_12*max_lambda_11 + CF_12)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11 + CF_14) +
      ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - 2*(-CS_14*max_lambda_11
      + CF_14) + ((3.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 +
      CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) - (-CS_14*max_lambda_11 +
      CF_14))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) + ((1.0/2.0))*(-CS_15*max_lambda_11 + CF_15) -
      (-CS_14*max_lambda_11 + CF_14)));

    beta_3 = ((1.0/36.0))*((-(-CS_10*max_lambda_11 + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) +
      ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) + ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))*(-(-CS_10*max_lambda_11
      + CF_10) - 9*(-CS_12*max_lambda_11 + CF_12) + ((9.0/2.0))*(-CS_11*max_lambda_11 + CF_11) +
      ((11.0/2.0))*(-CS_13*max_lambda_11 + CF_13))) + ((13.0/12.0))*((2*(-CS_11*max_lambda_11 + CF_11) -
      (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) - (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 +
      CF_13)*(2*(-CS_11*max_lambda_11 + CF_11) - (5.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) - CS_13*max_lambda_11 + CF_13)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 +
      CF_11))*(((1.0/2.0))*(-CS_13*max_lambda_11 + CF_13) - (3.0/2.0)*(-CS_12*max_lambda_11 + CF_12) -
      (1.0/2.0)*(-CS_10*max_lambda_11 + CF_10) + ((3.0/2.0))*(-CS_11*max_lambda_11 + CF_11)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_1 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/6.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((5.0/12.0))*(-CS_12*max_lambda_11 + CF_12))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_11*max_lambda_11 + CF_11) + ((1.0/8.0))*(-CS_13*max_lambda_11 + CF_13) +
      ((1.0/24.0))*(-CS_10*max_lambda_11 + CF_10) + ((13.0/24.0))*(-CS_12*max_lambda_11 + CF_12))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_15*max_lambda_11 + CF_15) +
      ((11.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_14*max_lambda_11 + CF_14) + ((1.0/6.0))*(-CS_12*max_lambda_11 + CF_12) +
      ((5.0/12.0))*(-CS_13*max_lambda_11 + CF_13))*delta_0*inv_omega_sum + Recon_1;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 +
      CF_23))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) - (1.0/2.0)*(CS_23*max_lambda_22 + CF_23))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 + CF_23) -
      (CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(CS_23*max_lambda_22 +
      CF_23) - (CS_22*max_lambda_22 + CF_22)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 + CF_23) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - 2*(CS_23*max_lambda_22 +
      CF_23) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(CS_22*max_lambda_22 +
      CF_22) + ((1.0/2.0))*(CS_24*max_lambda_22 + CF_24) - (CS_23*max_lambda_22 + CF_23)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 + CF_21) +
      ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))*(((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) - 2*(CS_21*max_lambda_22 +
      CF_21) + ((3.0/2.0))*(CS_22*max_lambda_22 + CF_22))) + ((13.0/12.0))*((((1.0/2.0))*(CS_20*max_lambda_22 + CF_20) +
      ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21))*(((1.0/2.0))*(CS_20*max_lambda_22 +
      CF_20) + ((1.0/2.0))*(CS_22*max_lambda_22 + CF_22) - (CS_21*max_lambda_22 + CF_21)));

    beta_3 = -(781.0/480.0)*(CS_24*max_lambda_22 + CF_24) - (781.0/1440.0)*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/36.0))*((9*(CS_23*max_lambda_22 + CF_23) - (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) -
      (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25)*(9*(CS_23*max_lambda_22 + CF_23) -
      (11.0/2.0)*(CS_22*max_lambda_22 + CF_22) - (9.0/2.0)*(CS_24*max_lambda_22 + CF_24) + CS_25*max_lambda_22 + CF_25))
      + ((13.0/12.0))*((2*(CS_24*max_lambda_22 + CF_24) - (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) -
      (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22)*(2*(CS_24*max_lambda_22 + CF_24) -
      (5.0/2.0)*(CS_23*max_lambda_22 + CF_23) - (1.0/2.0)*(CS_25*max_lambda_22 + CF_25) + CS_22*max_lambda_22 + CF_22))
      + ((781.0/480.0))*(CS_23*max_lambda_22 + CF_23) + ((781.0/1440.0))*(CS_25*max_lambda_22 + CF_25);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(CS_23*max_lambda_22 + CF_23))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_24*max_lambda_22 + CF_24) + ((1.0/8.0))*(CS_22*max_lambda_22 + CF_22) +
      ((1.0/24.0))*(CS_25*max_lambda_22 + CF_25) + ((13.0/24.0))*(CS_23*max_lambda_22 + CF_23))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_20*max_lambda_22 + CF_20) +
      ((11.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(CS_22*max_lambda_22 + CF_22))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_24*max_lambda_22 + CF_24))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 + CF_24) -
      (-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_22*max_lambda_22 + CF_22) + ((1.0/2.0))*(-CS_24*max_lambda_22 +
      CF_24) - (-CS_23*max_lambda_22 + CF_23)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22 + CF_22) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) - 2*(-CS_22*max_lambda_22
      + CF_22) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (-CS_22*max_lambda_22 +
      CF_22))*(((1.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) -
      (-CS_22*max_lambda_22 + CF_22)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22 + CF_24) +
      ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - 2*(-CS_24*max_lambda_22
      + CF_24) + ((3.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 +
      CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) - (-CS_24*max_lambda_22 +
      CF_24))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) + ((1.0/2.0))*(-CS_25*max_lambda_22 + CF_25) -
      (-CS_24*max_lambda_22 + CF_24)));

    beta_3 = ((1.0/36.0))*((-(-CS_20*max_lambda_22 + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) +
      ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) + ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))*(-(-CS_20*max_lambda_22
      + CF_20) - 9*(-CS_22*max_lambda_22 + CF_22) + ((9.0/2.0))*(-CS_21*max_lambda_22 + CF_21) +
      ((11.0/2.0))*(-CS_23*max_lambda_22 + CF_23))) + ((13.0/12.0))*((2*(-CS_21*max_lambda_22 + CF_21) -
      (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) - (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 +
      CF_23)*(2*(-CS_21*max_lambda_22 + CF_21) - (5.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) - CS_23*max_lambda_22 + CF_23)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 +
      CF_21))*(((1.0/2.0))*(-CS_23*max_lambda_22 + CF_23) - (3.0/2.0)*(-CS_22*max_lambda_22 + CF_22) -
      (1.0/2.0)*(-CS_20*max_lambda_22 + CF_20) + ((3.0/2.0))*(-CS_21*max_lambda_22 + CF_21)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_2 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/6.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((5.0/12.0))*(-CS_22*max_lambda_22 + CF_22))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_21*max_lambda_22 + CF_21) + ((1.0/8.0))*(-CS_23*max_lambda_22 + CF_23) +
      ((1.0/24.0))*(-CS_20*max_lambda_22 + CF_20) + ((13.0/24.0))*(-CS_22*max_lambda_22 + CF_22))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_25*max_lambda_22 + CF_25) +
      ((11.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_24*max_lambda_22 + CF_24) + ((1.0/6.0))*(-CS_22*max_lambda_22 + CF_22) +
      ((5.0/12.0))*(-CS_23*max_lambda_22 + CF_23))*delta_0*inv_omega_sum + Recon_2;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 +
      CF_33))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) - (1.0/2.0)*(CS_33*max_lambda_33 + CF_33))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 + CF_33) -
      (CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(CS_33*max_lambda_33 +
      CF_33) - (CS_32*max_lambda_33 + CF_32)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 + CF_33) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - 2*(CS_33*max_lambda_33 +
      CF_33) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(CS_32*max_lambda_33 +
      CF_32) + ((1.0/2.0))*(CS_34*max_lambda_33 + CF_34) - (CS_33*max_lambda_33 + CF_33)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 + CF_31) +
      ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))*(((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) - 2*(CS_31*max_lambda_33 +
      CF_31) + ((3.0/2.0))*(CS_32*max_lambda_33 + CF_32))) + ((13.0/12.0))*((((1.0/2.0))*(CS_30*max_lambda_33 + CF_30) +
      ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31))*(((1.0/2.0))*(CS_30*max_lambda_33 +
      CF_30) + ((1.0/2.0))*(CS_32*max_lambda_33 + CF_32) - (CS_31*max_lambda_33 + CF_31)));

    beta_3 = -(781.0/480.0)*(CS_34*max_lambda_33 + CF_34) - (781.0/1440.0)*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/36.0))*((9*(CS_33*max_lambda_33 + CF_33) - (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) -
      (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35)*(9*(CS_33*max_lambda_33 + CF_33) -
      (11.0/2.0)*(CS_32*max_lambda_33 + CF_32) - (9.0/2.0)*(CS_34*max_lambda_33 + CF_34) + CS_35*max_lambda_33 + CF_35))
      + ((13.0/12.0))*((2*(CS_34*max_lambda_33 + CF_34) - (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) -
      (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32)*(2*(CS_34*max_lambda_33 + CF_34) -
      (5.0/2.0)*(CS_33*max_lambda_33 + CF_33) - (1.0/2.0)*(CS_35*max_lambda_33 + CF_35) + CS_32*max_lambda_33 + CF_32))
      + ((781.0/480.0))*(CS_33*max_lambda_33 + CF_33) + ((781.0/1440.0))*(CS_35*max_lambda_33 + CF_35);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(CS_33*max_lambda_33 + CF_33))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_34*max_lambda_33 + CF_34) + ((1.0/8.0))*(CS_32*max_lambda_33 + CF_32) +
      ((1.0/24.0))*(CS_35*max_lambda_33 + CF_35) + ((13.0/24.0))*(CS_33*max_lambda_33 + CF_33))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_30*max_lambda_33 + CF_30) +
      ((11.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(CS_32*max_lambda_33 + CF_32))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_34*max_lambda_33 + CF_34))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 + CF_34) -
      (-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_32*max_lambda_33 + CF_32) + ((1.0/2.0))*(-CS_34*max_lambda_33 +
      CF_34) - (-CS_33*max_lambda_33 + CF_33)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33 + CF_32) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) - 2*(-CS_32*max_lambda_33
      + CF_32) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (-CS_32*max_lambda_33 +
      CF_32))*(((1.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) -
      (-CS_32*max_lambda_33 + CF_32)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33 + CF_34) +
      ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - 2*(-CS_34*max_lambda_33
      + CF_34) + ((3.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 +
      CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) - (-CS_34*max_lambda_33 +
      CF_34))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) + ((1.0/2.0))*(-CS_35*max_lambda_33 + CF_35) -
      (-CS_34*max_lambda_33 + CF_34)));

    beta_3 = ((1.0/36.0))*((-(-CS_30*max_lambda_33 + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) +
      ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) + ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))*(-(-CS_30*max_lambda_33
      + CF_30) - 9*(-CS_32*max_lambda_33 + CF_32) + ((9.0/2.0))*(-CS_31*max_lambda_33 + CF_31) +
      ((11.0/2.0))*(-CS_33*max_lambda_33 + CF_33))) + ((13.0/12.0))*((2*(-CS_31*max_lambda_33 + CF_31) -
      (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) - (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 +
      CF_33)*(2*(-CS_31*max_lambda_33 + CF_31) - (5.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) - CS_33*max_lambda_33 + CF_33)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 +
      CF_31))*(((1.0/2.0))*(-CS_33*max_lambda_33 + CF_33) - (3.0/2.0)*(-CS_32*max_lambda_33 + CF_32) -
      (1.0/2.0)*(-CS_30*max_lambda_33 + CF_30) + ((3.0/2.0))*(-CS_31*max_lambda_33 + CF_31)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_3 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/6.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((5.0/12.0))*(-CS_32*max_lambda_33 + CF_32))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_31*max_lambda_33 + CF_31) + ((1.0/8.0))*(-CS_33*max_lambda_33 + CF_33) +
      ((1.0/24.0))*(-CS_30*max_lambda_33 + CF_30) + ((13.0/24.0))*(-CS_32*max_lambda_33 + CF_32))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_35*max_lambda_33 + CF_35) +
      ((11.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_34*max_lambda_33 + CF_34) + ((1.0/6.0))*(-CS_32*max_lambda_33 + CF_32) +
      ((5.0/12.0))*(-CS_33*max_lambda_33 + CF_33))*delta_0*inv_omega_sum + Recon_3;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 +
      CF_43))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) - (1.0/2.0)*(CS_43*max_lambda_44 + CF_43))) +
      ((13.0/12.0))*((((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 + CF_43) -
      (CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(CS_43*max_lambda_44 +
      CF_43) - (CS_42*max_lambda_44 + CF_42)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 + CF_43) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - 2*(CS_43*max_lambda_44 +
      CF_43) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(CS_42*max_lambda_44 +
      CF_42) + ((1.0/2.0))*(CS_44*max_lambda_44 + CF_44) - (CS_43*max_lambda_44 + CF_43)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 + CF_41) +
      ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))*(((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) - 2*(CS_41*max_lambda_44 +
      CF_41) + ((3.0/2.0))*(CS_42*max_lambda_44 + CF_42))) + ((13.0/12.0))*((((1.0/2.0))*(CS_40*max_lambda_44 + CF_40) +
      ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41))*(((1.0/2.0))*(CS_40*max_lambda_44 +
      CF_40) + ((1.0/2.0))*(CS_42*max_lambda_44 + CF_42) - (CS_41*max_lambda_44 + CF_41)));

    beta_3 = -(781.0/480.0)*(CS_44*max_lambda_44 + CF_44) - (781.0/1440.0)*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/36.0))*((9*(CS_43*max_lambda_44 + CF_43) - (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) -
      (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45)*(9*(CS_43*max_lambda_44 + CF_43) -
      (11.0/2.0)*(CS_42*max_lambda_44 + CF_42) - (9.0/2.0)*(CS_44*max_lambda_44 + CF_44) + CS_45*max_lambda_44 + CF_45))
      + ((13.0/12.0))*((2*(CS_44*max_lambda_44 + CF_44) - (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) -
      (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42)*(2*(CS_44*max_lambda_44 + CF_44) -
      (5.0/2.0)*(CS_43*max_lambda_44 + CF_43) - (1.0/2.0)*(CS_45*max_lambda_44 + CF_45) + CS_42*max_lambda_44 + CF_42))
      + ((781.0/480.0))*(CS_43*max_lambda_44 + CF_43) + ((781.0/1440.0))*(CS_45*max_lambda_44 + CF_45);

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(CS_43*max_lambda_44 + CF_43))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(CS_44*max_lambda_44 + CF_44) + ((1.0/8.0))*(CS_42*max_lambda_44 + CF_42) +
      ((1.0/24.0))*(CS_45*max_lambda_44 + CF_45) + ((13.0/24.0))*(CS_43*max_lambda_44 + CF_43))*delta_3*inv_omega_sum +
      ((27.0/500.0))*(-(7.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_40*max_lambda_44 + CF_40) +
      ((11.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(CS_42*max_lambda_44 + CF_42))*delta_0*inv_omega_sum + Recon_4;

    beta_0 = ((1.0/4.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_44*max_lambda_44 + CF_44))) +
      ((13.0/12.0))*((((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 + CF_44) -
      (-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_42*max_lambda_44 + CF_42) + ((1.0/2.0))*(-CS_44*max_lambda_44 +
      CF_44) - (-CS_43*max_lambda_44 + CF_43)));

    beta_1 = ((1.0/4.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44 + CF_42) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) - 2*(-CS_42*max_lambda_44
      + CF_42) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (-CS_42*max_lambda_44 +
      CF_42))*(((1.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) -
      (-CS_42*max_lambda_44 + CF_42)));

    beta_2 = ((1.0/4.0))*((((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44 + CF_44) +
      ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - 2*(-CS_44*max_lambda_44
      + CF_44) + ((3.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 +
      CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) - (-CS_44*max_lambda_44 +
      CF_44))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) + ((1.0/2.0))*(-CS_45*max_lambda_44 + CF_45) -
      (-CS_44*max_lambda_44 + CF_44)));

    beta_3 = ((1.0/36.0))*((-(-CS_40*max_lambda_44 + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) +
      ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) + ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))*(-(-CS_40*max_lambda_44
      + CF_40) - 9*(-CS_42*max_lambda_44 + CF_42) + ((9.0/2.0))*(-CS_41*max_lambda_44 + CF_41) +
      ((11.0/2.0))*(-CS_43*max_lambda_44 + CF_43))) + ((13.0/12.0))*((2*(-CS_41*max_lambda_44 + CF_41) -
      (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) - (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 +
      CF_43)*(2*(-CS_41*max_lambda_44 + CF_41) - (5.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) - CS_43*max_lambda_44 + CF_43)) +
      ((781.0/720.0))*((((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 +
      CF_41))*(((1.0/2.0))*(-CS_43*max_lambda_44 + CF_43) - (3.0/2.0)*(-CS_42*max_lambda_44 + CF_42) -
      (1.0/2.0)*(-CS_40*max_lambda_44 + CF_40) + ((3.0/2.0))*(-CS_41*max_lambda_44 + CF_41)));

   inv_beta_0 = 1.0/(eps + beta_0);

   inv_beta_1 = 1.0/(eps + beta_1);

   inv_beta_2 = 1.0/(eps + beta_2);

   inv_beta_3 = 1.0/(eps + beta_3);

    alpha_0 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_0));

    alpha_1 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_1));

    alpha_2 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_2));

    alpha_3 = ((1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 +
      fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 -
      (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 +
      ((1.0/6.0))*beta_0 + ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3)*(1 + fabs(-beta_3 - (2.0/3.0)*beta_1 + ((1.0/6.0))*beta_0 +
      ((1.0/6.0))*beta_2)*inv_beta_3));

   inv_alpha_sum = 1.0/((alpha_0 + alpha_1 + alpha_2 + alpha_3));

   delta_0 = ((alpha_0*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_1 = ((alpha_1*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_2 = ((alpha_2*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

   delta_3 = ((alpha_3*inv_alpha_sum < TENO_CT) ? (
   0.0
)
: (
   1.0
));

    inv_omega_sum = 1.0/((((3.0/10.0))*delta_1 + ((23.0/125.0))*delta_3 + ((27.0/500.0))*delta_2 +
      ((231.0/500.0))*delta_0));

    Recon_4 = ((3.0/10.0))*(-(1.0/12.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/6.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((5.0/12.0))*(-CS_42*max_lambda_44 + CF_42))*delta_1*inv_omega_sum +
      ((23.0/125.0))*(-(5.0/24.0)*(-CS_41*max_lambda_44 + CF_41) + ((1.0/8.0))*(-CS_43*max_lambda_44 + CF_43) +
      ((1.0/24.0))*(-CS_40*max_lambda_44 + CF_40) + ((13.0/24.0))*(-CS_42*max_lambda_44 + CF_42))*delta_3*inv_omega_sum
      + ((27.0/500.0))*(-(7.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_45*max_lambda_44 + CF_45) +
      ((11.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_2*inv_omega_sum +
      ((231.0/500.0))*(-(1.0/12.0)*(-CS_44*max_lambda_44 + CF_44) + ((1.0/6.0))*(-CS_42*max_lambda_44 + CF_42) +
      ((5.0/12.0))*(-CS_43*max_lambda_44 + CF_43))*delta_0*inv_omega_sum + Recon_4;

    wk10_B0(0,0,0) = 0.707106781186547*AVG_2_rho*Recon_3*inv_AVG_a + 0.707106781186547*AVG_2_rho*Recon_4*inv_AVG_a +
      Recon_2;

    wk11_B0(0,0,0) = AVG_2_u0*Recon_2 - AVG_2_rho*Recon_1 + 0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_3*inv_AVG_a +
      0.707106781186547*AVG_2_rho*AVG_2_u0*Recon_4*inv_AVG_a;

    wk12_B0(0,0,0) = AVG_2_rho*Recon_0 + AVG_2_u1*Recon_2 + 0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_3*inv_AVG_a +
      0.707106781186547*AVG_2_rho*AVG_2_u1*Recon_4*inv_AVG_a;

    wk13_B0(0,0,0) = AVG_2_u2*Recon_2 + 0.707106781186547*(-AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a +
      0.707106781186547*(AVG_2_a + AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a;

    wk14_B0(0,0,0) = (((1.0/2.0))*(AVG_2_u0*AVG_2_u0) + ((1.0/2.0))*(AVG_2_u1*AVG_2_u1) +
      ((1.0/2.0))*(AVG_2_u2*AVG_2_u2))*Recon_2 + AVG_2_rho*AVG_2_u1*Recon_0 - AVG_2_rho*AVG_2_u0*Recon_1 +
      0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
      (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 + AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_3*inv_AVG_a +
      0.707106781186547*(((AVG_2_a*AVG_2_a) + ((1.0/2.0))*((AVG_2_u0*AVG_2_u0) + (AVG_2_u1*AVG_2_u1) +
      (AVG_2_u2*AVG_2_u2))*gamma_m1)*invgamma_m1 - AVG_2_a*AVG_2_u2)*AVG_2_rho*Recon_4*inv_AVG_a;

}

 void opensbliblock00Kernel013(const ACC<double> &wk0_B0, const ACC<double> &wk10_B0, const ACC<double> &wk11_B0, const
ACC<double> &wk12_B0, const ACC<double> &wk13_B0, const ACC<double> &wk14_B0, const ACC<double> &wk1_B0, const
ACC<double> &wk2_B0, const ACC<double> &wk3_B0, const ACC<double> &wk4_B0, const ACC<double> &wk5_B0, const ACC<double>
&wk6_B0, const ACC<double> &wk7_B0, const ACC<double> &wk8_B0, const ACC<double> &wk9_B0, ACC<double> &Residual0_B0,
ACC<double> &Residual1_B0, ACC<double> &Residual2_B0, ACC<double> &Residual3_B0, ACC<double> &Residual4_B0)
{
    Residual0_B0(0,0,0) = -(-wk0_B0(-1,0,0) + wk0_B0(0,0,0))*invDelta0block0 - (-wk5_B0(0,-1,0) +
      wk5_B0(0,0,0))*invDelta1block0 - (-wk10_B0(0,0,-1) + wk10_B0(0,0,0))*invDelta2block0;

    Residual1_B0(0,0,0) = -(-wk1_B0(-1,0,0) + wk1_B0(0,0,0))*invDelta0block0 - (-wk6_B0(0,-1,0) +
      wk6_B0(0,0,0))*invDelta1block0 - (-wk11_B0(0,0,-1) + wk11_B0(0,0,0))*invDelta2block0;

    Residual2_B0(0,0,0) = -(-wk2_B0(-1,0,0) + wk2_B0(0,0,0))*invDelta0block0 - (-wk7_B0(0,-1,0) +
      wk7_B0(0,0,0))*invDelta1block0 - (-wk12_B0(0,0,-1) + wk12_B0(0,0,0))*invDelta2block0;

    Residual3_B0(0,0,0) = -(-wk3_B0(-1,0,0) + wk3_B0(0,0,0))*invDelta0block0 - (-wk8_B0(0,-1,0) +
      wk8_B0(0,0,0))*invDelta1block0 - (-wk13_B0(0,0,-1) + wk13_B0(0,0,0))*invDelta2block0;

    Residual4_B0(0,0,0) = -(-wk4_B0(-1,0,0) + wk4_B0(0,0,0))*invDelta0block0 - (-wk9_B0(0,-1,0) +
      wk9_B0(0,0,0))*invDelta1block0 - (-wk14_B0(0,0,-1) + wk14_B0(0,0,0))*invDelta2block0;

}

void opensbliblock00Kernel016(const ACC<double> &u0_B0, ACC<double> &wk0_B0)
{
    wk0_B0(0,0,0) = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel018(const ACC<double> &u1_B0, ACC<double> &wk1_B0)
{
    wk1_B0(0,0,0) = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel020(const ACC<double> &u2_B0, ACC<double> &wk2_B0)
{
    wk2_B0(0,0,0) = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

}

void opensbliblock00Kernel022(const ACC<double> &u0_B0, ACC<double> &wk3_B0)
{
    wk3_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
      ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel023(const ACC<double> &u1_B0, ACC<double> &wk4_B0)
{
    wk4_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
      ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel024(const ACC<double> &u2_B0, ACC<double> &wk5_B0)
{
    wk5_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
      ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

}

void opensbliblock00Kernel025(const ACC<double> &u0_B0, ACC<double> &wk6_B0)
{
    wk6_B0(0,0,0) = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel026(const ACC<double> &u1_B0, ACC<double> &wk7_B0)
{
    wk7_B0(0,0,0) = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

}

void opensbliblock00Kernel027(const ACC<double> &u2_B0, ACC<double> &wk8_B0)
{
    wk8_B0(0,0,0) = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

}

 void opensbliblock00Kernel035(const ACC<double> &T_B0, const ACC<double> &mu_B0, const ACC<double> &u0_B0, const
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

 void opensbliblock00Kernel049(const ACC<double> &Residual0_B0, const ACC<double> &Residual1_B0, const ACC<double>
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

 void opensbliblock00Kernel048(const ACC<double> &mu_B0, const ACC<double> &rho_B0, const ACC<double> &u0_B0, const
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
    d1_u1_dz = (-(2.0/3.0)*u1_B0(0,0,-1) - (1.0/12.0)*u1_B0(0,0,2) + ((1.0/12.0))*u1_B0(0,0,-2) +
      ((2.0/3.0))*u1_B0(0,0,1))*invDelta2block0;

    d1_u2_dy = (-(2.0/3.0)*u2_B0(0,-1,0) - (1.0/12.0)*u2_B0(0,2,0) + ((1.0/12.0))*u2_B0(0,-2,0) +
      ((2.0/3.0))*u2_B0(0,1,0))*invDelta1block0;

   wx = -d1_u1_dz + d1_u2_dy;

    d1_u2_dx = (-(2.0/3.0)*u2_B0(-1,0,0) - (1.0/12.0)*u2_B0(2,0,0) + ((1.0/12.0))*u2_B0(-2,0,0) +
      ((2.0/3.0))*u2_B0(1,0,0))*invDelta0block0;

    d1_u0_dz = (-(2.0/3.0)*u0_B0(0,0,-1) - (1.0/12.0)*u0_B0(0,0,2) + ((1.0/12.0))*u0_B0(0,0,-2) +
      ((2.0/3.0))*u0_B0(0,0,1))*invDelta2block0;

   wy = -d1_u2_dx + d1_u0_dz;

    d1_u1_dx = (-(2.0/3.0)*u1_B0(-1,0,0) - (1.0/12.0)*u1_B0(2,0,0) + ((1.0/12.0))*u1_B0(-2,0,0) +
      ((2.0/3.0))*u1_B0(1,0,0))*invDelta0block0;

    d1_u0_dy = (-(2.0/3.0)*u0_B0(0,-1,0) - (1.0/12.0)*u0_B0(0,2,0) + ((1.0/12.0))*u0_B0(0,-2,0) +
      ((2.0/3.0))*u0_B0(0,1,0))*invDelta1block0;

   wz = -d1_u0_dy + d1_u1_dx;

    d1_u0_dx = (-(2.0/3.0)*u0_B0(-1,0,0) - (1.0/12.0)*u0_B0(2,0,0) + ((1.0/12.0))*u0_B0(-2,0,0) +
      ((2.0/3.0))*u0_B0(1,0,0))*invDelta0block0;

    d1_u2_dz = (-(2.0/3.0)*u2_B0(0,0,-1) - (1.0/12.0)*u2_B0(0,0,2) + ((1.0/12.0))*u2_B0(0,0,-2) +
      ((2.0/3.0))*u2_B0(0,0,1))*invDelta2block0;

    d1_u1_dy = (-(2.0/3.0)*u1_B0(0,-1,0) - (1.0/12.0)*u1_B0(0,2,0) + ((1.0/12.0))*u1_B0(0,-2,0) +
      ((2.0/3.0))*u1_B0(0,1,0))*invDelta1block0;

   divV_B0(0,0,0) = d1_u0_dx + d1_u1_dy + d1_u2_dz;

   *rhom_B0 = rho_B0(0,0,0) + *rhom_B0;

    *KE_B0 = 0.5*((u0_B0(0,0,0)*u0_B0(0,0,0)) + (u1_B0(0,0,0)*u1_B0(0,0,0)) + (u2_B0(0,0,0)*u2_B0(0,0,0)))*rho_B0(0,0,0)
      + *KE_B0;

   *dilatation_dissipation_B0 = ((4.0/3.0))*(divV_B0(0,0,0)*divV_B0(0,0,0))*mu_B0(0,0,0) + *dilatation_dissipation_B0;

   *enstrophy_dissipation_B0 = ((wx*wx) + (wy*wy) + (wz*wz))*mu_B0(0,0,0) + *enstrophy_dissipation_B0;

}

#endif
