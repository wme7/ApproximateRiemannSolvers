function [qr,ql] = weno5_char(q,maxnx,num_eqn,num_ghost,evl,evr)

% This one uses characteristic decomposition
% evl, evr are left and right eigenvectors at each interface

% integer,          intent(in) :: maxnx, num_eqn, num_ghost
% double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
% double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost)
% double precision, intent(out) :: qr(num_eqn,maxnx+2*num_ghost)
% double precision, intent(in) :: evl(num_eqn,num_eqn,maxnx+2*num_ghost)
% double precision, intent(in) :: evr(num_eqn,num_eqn,maxnx+2*num_ghost)

% integer :: mx2

mx2 = size(q,2);

% loop over all equations (all components).  
% the reconstruction is performed using characteristic decomposition

for m=1:num_eqn
    for i=2:mx2
        % compute and store the differences of the cell averages
        dq(m,i)=q(m,i)-q(m,i-1);
    end
end 

for m=1:num_eqn
    for i=num_ghost:mx2-num_ghost+1
        % Compute the part of the reconstruction that is stencil-independent
        qr(m,i-1) = (-q(m,i-2)+7.*(q(m,i-1)+q(m,i))-q(m,i+1))/12.;
        ql(m,i)   = qr(m,i-1);
    end
end 

for ip=1:num_eqn
    
    % Project the difference of the cell averages to the 'm'th characteristic field
    for m2 = -2:2
       for  i = num_ghost+1:mx2-2
          hh(m2,i) = 0.d0;
          for m=1:num_eqn 
            hh(m2,i) = hh(m2,i) + evl(ip,m,i)*dq(m,i+m2);
          end
       end
    end

    % the reconstruction
    for m1=1:2

        % m1=1: construct ql
        % m1=2: construct qr

        im=(-1)^(m1+1);
        ione=im;
        inone=-im;
        intwo=-2*im;

        for i=num_ghost:mx2-num_ghost+1

            t1=im*(hh(intwo,i)-hh(inone,i));
            t2=im*(hh(inone,i)-hh(0,i    ));
            t3=im*(hh(0,i    )-hh(ione,i ));

            tt1=13.*t1^2+3.*(   hh(intwo,i)-3.*hh(inone,i))^2;
            tt2=13.*t2^2+3.*(   hh(inone,i)+   hh(0,i    ))^2;
            tt3=13.*t3^2+3.*(3.*hh(0,i    )-   hh(ione,i ))^2;

            tt1=(epweno+tt1)^2;
            tt2=(epweno+tt2)^2;
            tt3=(epweno+tt3)^2;
            s1 =tt2*tt3;
            s2 =6.*tt1*tt3;
            s3 =3.*tt1*tt2;
            t0 =1./(s1+s2+s3);
            s1 =s1*t0;
            s3 =s3*t0;

            uu(m1,i) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.;

        end %end loop over interfaces
    end %end loop over which side of interface

    % Project to the physical space:
    for m = 1:num_eqn
        for i=num_ghost:mx2-num_ghost+1
            qr(m,i-1) = qr(m,i-1) + evr(m,ip,i)*uu(1,i);
            ql(m,i  ) = ql(m,i  ) + evr(m,ip,i)*uu(2,i);
        end 
    end
end 

end %subroutine weno5_char