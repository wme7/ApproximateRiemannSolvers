function [ql,qr]= weno5(q)

% Traditional WENO5 reconstruction

% double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
% double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)
qr=zeros(size(q));
ql=zeros(size(q));
epweno=1E-40;

% integer :: num_eqn, mx2

num_ghost = 3;
num_eqn = size(q,1);
mx2 = size(q,2); 

% loop over all equations (all components).  
% the reconstruction is performed component-wise;
% no characteristic decomposition is used here

for m=1:num_eqn

    for i=2:mx2
        % compute and store the differences of the cell averages
        dq1m(i) = q(m,i)-q(m,i-1);
    end 

    % the reconstruction

    for m1=1:2

        % m1=1: construct ql
        % m1=2: construct qr

        im=(-1)^(m1+1);
        ione=im; inone=-im; intwo=-2*im;

        for i=num_ghost:mx2-num_ghost+1

            t1=im*(dq1m(i+intwo)-dq1m(i+inone));
            t2=im*(dq1m(i+inone)-dq1m(i      ));
            t3=im*(dq1m(i      )-dq1m(i+ione ));

            tt1=13.*t1^2+3.*(   dq1m(i+intwo)-3.*dq1m(i+inone))^2;
            tt2=13.*t2^2+3.*(   dq1m(i+inone)+   dq1m(i      ))^2;
            tt3=13.*t3^2+3.*(3.*dq1m(i      )-   dq1m(i+ione ))^2;

            tt1=(epweno+tt1)^2;
            tt2=(epweno+tt2)^2;
            tt3=(epweno+tt3)^2;
            s1 =tt2*tt3;
            s2 =6.*tt1*tt3;
            s3 =3.*tt1*tt2;
            t0 =1./(s1+s2+s3);
            s1 =s1*t0;
            s3 =s3*t0;

            uu(m1,i) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3 ...
                     +(-q(m,i-2)+7*(q(m,i-1)+q(m,i))-q(m,i+1))/12;
        end 
    end 

   qr(m,num_ghost-1:mx2-num_ghost  )=uu(1,num_ghost:mx2-num_ghost+1);
   ql(m,num_ghost  :mx2-num_ghost+1)=uu(2,num_ghost:mx2-num_ghost+1);

end %subroutine weno5