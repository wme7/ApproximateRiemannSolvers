%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 2-D Euler system of equations with 5th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%       dq_i/dt + df_i/dx = 0, for x \in [a,b]^2 and i =1,...,D
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2012.08.25
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coded by Manuel A. Diaz, 2012.12.27. Last modif: 29.04.2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: C.-W. Shu, High order weighted essentially non-oscillatory schemes
% for convection dominated problems, SIAM Review, 51:82-126, (2009). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 
% 1. A fully conservative finite volume implementation of the method of
% lines (MOL) using WENO5 associated with SSP-RK33 time integration method. 
% 2. Sharpenning of contact discontinuities is NOT implemented here.

clear; %close all; clc;
global gamma

%% Parameters
CFL     = 0.55;	% CFL number
tFinal	= 0.10;	% Final time
nx      = 200;  % Number of cells
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = 01;	% 10 IC cases are available
plot_fig= 1;

% Discretize spatial domain
Lx=1; dx=Lx/nx; xc=dx/2:dx:Lx;

% Set IC
[r0,u0,p0] = Euler_IC1d(xc,IC);
E0 = p0./((gamma-1)*r0)+0.5*u0.^2;  % Total Energy density
a0 = sqrt(gamma*p0./r0);            % Speed of sound
q0=[r0; r0.*u0; r0.*E0];        % vec. of conserved properties

% Exact solution
[xe,re,ue,pe,ee,te,Me,se] = ...
   EulerExact(r0(1),u0(1),p0(1),r0(nx),u0(nx),p0(nx),tFinal,gamma);
Ee = pe./((gamma-1)*re)+0.5*ue.^2;

% Adjust grid for ghost cells
nx=nx+4; zero=[0;0;0]; q0=[zero,zero,q0,zero,zero];

% Boundary Conditions in ghost cells
q0(:,1 )=q0(:, 3  ); q0(:, 2  )=q0(:, 3  ); % Natural BCs
q0(:,nx)=q0(:,nx-2); q0(:,nx-1)=q0(:,nx-2); 

% Initial time step
lambda0=abs(u0)+a0; dt0=CFL*dx/max(lambda0(:));

% Load IC
q=q0; t=0; it=0; dt=dt0; lambda=lambda0;

%% Solver Loop
while t<tFinal
    % RK Initial step
    qo = q;
    
    % 1st stage
    dF=FV_WENO5LF_1d(q,max(lambda(:)),nx,dx);     q = qo-dt*dF;
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
    
    % 2nd Stage
    dF=FV_WENO5LF_1d(q,max(lambda(:)),nx,dx);     q = 0.75*qo+0.25*(q-dt*dF);
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
    
    % 3rd stage
    dF=FV_WENO5LF_1d(q,max(lambda(:)),nx,dx);     q = (qo+2*(q-dt*dF))/3;
    q(:,1)=qo(:,3); q(:, nx )=qo(:,nx-2); % Neumann BCs
    q(:,2)=qo(:,3); q(:,nx-1)=qo(:,nx-2); % Neumann BCs
    
    % compute flow properties
    r=q(1,:); u=q(2,:)./r; E=q(3,:)./r; p=(gamma-1)*r.*(E-0.5*u.^2); a=sqrt(gamma*p./r);
    
    % Update dt and time
    lambda=abs(u)+a; dt=CFL*dx/max(lambda(:));
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt; it=it+1;
    
    % Plot figure
    if rem(it,10) == 0
        if plot_fig == 1
            subplot(2,2,1); plot(xc,r(3:nx-2),'.b',xe,re);
            subplot(2,2,2); plot(xc,u(3:nx-2),'.m',xe,ue);
            subplot(2,2,3); plot(xc,p(3:nx-2),'.k',xe,pe);
            subplot(2,2,4); plot(xc,E(3:nx-2),'.r',xe,Ee);
            drawnow
        end
    end
end

% Remove ghost cells
q=q(:,3:nx-2); nx=nx-4; 

% compute flow properties
r=q(1,:); u=q(2,:)./r; E=q(3,:)./r; p=(gamma-1)*r.*(E-0.5*u.^2);

% Plots results
figure(1);
subplot(2,2,1); plot(xc,r,'ro',xe,re,'-k'); xlabel('x'); ylabel('\rho'); legend('WENO-LF','Exact'); 
title('SSP-RK3 FV-WENO-LF Euler Eqns.')
subplot(2,2,2); plot(xc,u,'ro',xe,ue,'-k'); xlabel('x'); ylabel('u'); 
subplot(2,2,3); plot(xc,p,'ro',xe,pe,'-k'); xlabel('x'); ylabel('p'); 
subplot(2,2,4); plot(xc,E,'ro',xe,Ee,'-k'); xlabel('x'); ylabel('E'); 