% example of reading a 2D file into matlab
% compressible laminar channel comparison with exact solution (constant
% viscosity)
clear all
file='opensbli_output.h5';

Delta0block0=h5read(file,'/Delta0block0'); 
Delta1block0=h5read(file,'/Delta1block0');
Minf=h5read(file,'/Minf');
Pr=h5read(file,'/Pr');
Re=h5read(file,'/Re');
Twall=h5read(file,'/Twall');
N0=h5read(file,'/block0np0');
N1=h5read(file,'/block0np1');
gamma=h5read(file,'/gama');
dt=h5read(file,'/dt');

dx=Delta0block0;
Lx=Delta0block0*double(N0);

fprintf('\nMinf, Re, Pr, Tw = %f %f %f %f\n',Minf,Re,Pr,Twall)
fprintf('Nx, Ny = %i %i \n',N0,N1)
time=h5read(file,'/simulation_time');
fprintf('\ntime = %f \n',time)

% read into dummy arrays and remove halos (keep the periodic element in x)
dum=h5read(file,'/opensbliblock00/x0_B0'); x=dum(6:N0+6,6:N1+5);
dum=h5read(file,'/opensbliblock00/x1_B0'); y=dum(6:N0+6,6:N1+5);
dum=h5read(char(file),'/opensbliblock00/rho_B0'); rho=dum(6:N0+6,6:N1+5);
dum=h5read(char(file),'/opensbliblock00/rhou0_B0'); rhou=dum(6:N0+6,6:N1+5);
dum=h5read(char(file),'/opensbliblock00/rhou1_B0'); rhov=dum(6:N0+6,6:N1+5);
dum=h5read(char(file),'/opensbliblock00/rhoE_B0'); rhoE=dum(6:N0+6,6:N1+5);

% deduced variables
u=rhou./rho;
v=rhov./rho;
e=rhoE./rho-0.5*(u.^2+v.^2);
p=(gamma-1.0)*rho.*e;
T=e.*Minf.^2*gamma*(gamma-1.0);
M=sqrt((u.^2+v.^2).*rho./(gamma.*p));

% some useful diagnositcs
xmin=min(min(min(x)));
xmax=max(max(max(x)));
ymin=min(min(min(y)));
ymax=max(max(max(y)));
rhomin=min(min(min(rho)));
rhomax=max(max(max(rho)));
umin=min(min(min(u)));
umax=max(max(max(u)));
vmin=min(min(min(v)));
vmax=max(max(max(v)));
pmin=min(min(min(p)));
pmax=max(max(max(p)));
Tmin=min(min(min(T)));
Tmax=max(max(max(T)));
Mmin=min(min(min(M)));
Mmax=max(max(max(M)));
fprintf('x min, x max = %f %f \n',xmin,xmax)
fprintf('y min, y max = %f %f \n',ymin,ymax)
fprintf('rho min, rho max = %f %f \n',rhomin,rhomax)
fprintf('u min, u max = %f %f \n',umin,umax)
fprintf('v min, v max = %f %f \n',vmin,vmax)
fprintf('p min, p max = %f %f \n',pmin,pmax)
fprintf('T min, T max = %f %f \n',Tmin,Tmax)
fprintf('M min, M max = %f %f \n',Mmin,Mmax)

% integrals of rho and rhou over whole domain
rsum=0.0;
msum=0.0;
for i=1:N0
    for j=1:N1-1
        rsum=rsum+dx*0.5*(y(i,j+1)-y(i,j))*(rho(i,j)+rho(i,j+1));
        msum=msum+dx*0.5*(y(i,j+1)-y(i,j))*(rho(i,j)*u(i,j)+rho(i,j+1)*u(i,j+1));
    end
end
rsum=rsum/Lx;
msum=msum/Lx;
ucl=u(1,(N1-1)/2+1,1);

% theoretical solution on actual grid
for j=1:N1
    ulam(j)=0.5*Re*(1.0-y(1,j)^2);
    Tlam(j)=1.0+(gamma-1.0)*Pr*Re^2*Minf^2*(1.0-y(1,j)^4)/12.0;
    rholam(j)=1.0/Tlam(j);
    plam(j)=1.0/(gamma*Minf^2);
end

% integrals of exact solution (on the actual grid)
rsume=0.0;
msume=0.0;
for j=1:N1-1
    rsume=rsume+0.5*(y(1,j+1,1)-y(1,j,1))*(rholam(j)+rholam(j+1));
    msume=msume+0.5*(y(1,j+1,1)-y(1,j,1))*(rholam(j)*ulam(j)+rholam(j+1)*ulam(j+1));
end

% centreline velocities
ucl=u(1,(N1-1)/2+1,1);
ulamcl=ulam(1,(N1-1)/2+1,1);

% print the results
fprintf('\nrelative error in centreline velocity = %g \n',(ucl-ulamcl)/ulamcl)
fprintf('relative error in integral of density = %g \n',(rsum-rsume)/rsume)
fprintf('relative error in mass flowrate = %g \n',(msum-msume)/msume)

% plot the actual vs exact u,rho,T
subplot(3,1,1)
plot(y(1,:),ulam,y(1,:),u(1,:),'o')
legend('Exact','OpenSBLI')
xlabel('y')
ylabel('u')

subplot(3,1,2)
plot(y(1,:),Tlam,y(1,:),T(1,:),'o')
legend('Exact','OpenSBLI')
xlabel('y')
ylabel('T')

subplot(3,1,3)
plot(y(1,:),rholam,y(1,:),rho(1,:),'o')
legend('Exact','OpenSBLI')
xlabel('y')
ylabel('\rho')

