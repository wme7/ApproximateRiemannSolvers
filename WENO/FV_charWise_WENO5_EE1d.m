function res = FV_charWise_WENO5_EE1d(a,q,~,dx)
% *************************************************************************
%
%    Characteristic-wise Finite Volume-1d for the Euler Equations 
%
% Based on:
% ---------
% [1] Shu, Chi-Wang. "Essentially non-oscillatory and weighted essentially 
%     non-oscillatory schemes for hyperbolic conservation laws." Advanced 
%     numerical approximation of nonlinear hyperbolic equations. Springer, 
%     Berlin, Heidelberg, 1998. 325-432.
% [2] Jiang, Guang-Shan, and Cheng-chin Wu. "A high-order WENO finite
%     difference scheme for the equations of ideal magnetohydrodynamics."
%     Journal of Computational Physics 150.2 (1999): 561-594.  
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
%       last updated on 2018.06.20, NHRI Taiwan.
% *************************************************************************
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |___________S2__________|
%                               |                       |
%                       |___________S1__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: the reconstruction is performed using characteristic decomposition
% NOTE: Roe averages are assumed for the properties at the cell interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gamma
E = size(q,1); % number of equations
N = size(q,2); % number of nodes
R=3; I=R:(N-R); % R: stencil size 
epweno=1E-40;
gamma1 = gamma-1;

qr=zeros(size(q));
ql=zeros(size(q));
LF=zeros(size(q)); 
res=zeros(size(q));

evr = zeros(3,3,N);
evl = zeros(3,3,N);

% Compute eigenvectors at the cell interfaces 
for i = 2:N
    % Compute properties at cell interfaces using Roe avegares
    r_sqrtl = sqrt(q(1,i-1));
    r_sqrtr = sqrt(q(1, i ));
    pl = gamma1*(q(3,i-1) - 0.5*(q(2,i-1)^2)/q(1,i-1));
    pr = gamma1*(q(3, i ) - 0.5*(q(2, i )^2)/q(1, i ));
    r_sq2 = r_sqrtl + r_sqrtr;
    u = (q(2,i-1)/r_sqrtl + q(2,i)/r_sqrtr)/r_sq2;
    H = (((q(3,i-1)+pl)/r_sqrtl + (q(3,i)+pr)/r_sqrtr))/r_sq2;
    c2 = gamma1*(H - 0.5*u^2);
    c = sqrt(c2);

    % Construct matrix of right eigenvectors
    %      _                    _ 
    %     |                      |
    %     |   1      1       1   |
    %     |                      |
    % R = |  u-c     u      u+c  |
    %     |                      |
    %     |  H-uc   u^2/2   H+uc |
    %     |_                    _|
    
    evr(:,:,i) = [...
          1  ,  1  ,  1   ;...
         u-c ,  u  , u+c  ;...
        H-u*c,u^2/2,H+u*c];

    % Construct matrix of left eigenvectors
    %                          _                                       _ 
    %                         |                                         |
    %                         |  uc/(gamma-1)+u^2/2  -c/(gamma-1)-u   1 |
    %                         |                                         |
    % R^{-1}=(gamma-1)/(2c^2)*|  2(H-u^2)             2u             -2 |
    %                         |                                         |
    %                         | -uc/(gamma-1)+u^2/2   c/(gamma-1)-u   1 |
    %                         |_                                       _|

    evl(:,:,i) = gamma1/(2*c^2)*[...
         c*u/gamma1+u^2/2,-(c/gamma1+u), 1 ;...
              2*(H-u^2)  ,    2*u      ,-2 ;...
        -c*u/gamma1+u^2/2, c/gamma1-u  , 1];
end 

% compute and store the differences of the cell averages
for i=2:N
    dqmh(:,i)=q(:,i)-q(:,i-1); % dq_{i-1/2}
end
    
% Compute the part of the reconstruction that is stencil-independent
for i=R:N-R+1
    qr(:,i-1) = (-q(:,i-2)+7.*(q(:,i-1)+q(:,i))-q(:,i+1))/12.;
    ql(:, i ) = qr(:,i-1);
end

% Produce the WENO reconstruction
for ip=1:E
    
    % Project the difference of the cell averages to the 'm'th
    % characteristic field: qs
    for m2 = -2:2
       for  i = R+1:N-2
          qs(m2+3,i) = 0;
          for e=1:E 
            qs(m2+3,i) = qs(m2+3,i) + evl(ip,e,i)*dqmh(e,i+m2);
          end
       end
    end

    % the reconstruction
    for idx=1:2

        % idx=1: construct hn (qr)
        % idx=2: construct hp (ql)

        im=(-1)^(idx+1);
        i1=im+R; in1=-im+R; in2=-2*im+R;

        for i=R:N-R+1

            t1=im*(qs(in2,i)-qs(in1,i));
            t2=im*(qs(in1,i)-qs(R, i ));
            t3=im*(qs(R, i )-qs(i1,i ));

            IS1=13.*t1^2+3.*(   qs(in2,i)-3.*qs(in1,i))^2;
            IS2=13.*t2^2+3.*(   qs(in1,i)+   qs(R, i ))^2;
            IS3=13.*t3^2+3.*(3.*qs(R, i )-   qs(i1,i ))^2;

            IS1=(epweno+IS1)^2;
            IS2=(epweno+IS2)^2;
            IS3=(epweno+IS3)^2;
            s1 =IS2*IS3;
            s2 =6.*IS1*IS3;
            s3 =3.*IS1*IS2;
            t0 =1./(s1+s2+s3);
            s1 =s1*t0;
            s3 =s3*t0;

            h(idx,i) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.;

        end % loop over interfaces
    end % loop over which side of interface

    % Project to the physical space:
    for e = 1:E
        for i=R:N-R+1
            qr(e,i-1) = qr(e,i-1) + evr(e,ip,i)*h(1,i);
            ql(e, i ) = ql(e, i ) + evr(e,ip,i)*h(2,i);
        end 
    end
    
end 

%% Compute finite volume residual term, df/dx.
LF(:,I) = 0.5*(F(qr(:,I))+F(ql(:,I+1))-abs(a).*(ql(:,I+1)-qr(:,I))); % Lax friedrichs flux
% for j = I % for all faces of the domain cells
%     res(:, j ) = res(:, j ) + LF(:,j)/dx;
%     res(:,j+1) = res(:,j+1) - LF(:,j)/dx;
% end % or alternatively :
res(:,I) = (LF(:,I)-LF(:,I-1))/dx; % L = -df(q)/dx.

% Flux contribution of the LEFT MOST FACE: left face of cell j=3.
res(:,3) = res(:,3)-LF(:,3)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-2.
res(:,N-2)=res(:,N-2)+LF(:,N-2)/dx;
end

% Compute flux vector
function flux = F(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:); p=(gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector of conserved properties
    flux=[rho.*u; rho.*u.^2+p; u.*(E+p)];
end