function R = eigenmatrix(rho,rhou,rhov,rhow,E,B1,B2,B3,n1,n2,n3)
R = zeros(8,8);
global gamma

u1 = rhou/rho;
u2 = rhov/rho;
u3 = rhow/rho;

p = (gamma - 1)*(E - 0.5*rho*(u1^2 + u2^2 + u3^2) - 0.5*(B1^2 + B2^2 + B3^2));


[t1,t2,t3] = compute_t_vector(n1, n2, n3, B1, B2, B3);

rhosq = sqrt(abs(rho));
g1=gamma-1;
a2 = abs(gamma*p/rho);
a = sqrt(a2);

sqg2=sqrt(1/(2*gamma));

sq12 = sqrt(0.5);

sqpr = sqrt(abs(p))/rho;

sqpor = sqrt(abs(p/rho));

sq1og = sqrt(1/gamma);

b1s=B1/rhosq;
b2s=B2/rhosq;
b3s=B3/rhosq;

bN = b1s*n1+b2s*n2+b3s*n3;

d = a2 + (b1s^2 + b2s^2 + b3s^2);
cf = sqrt(0.5d0*abs(d+sqrt(abs(d^2-4*a2*(bN)^2))));
cs = sqrt(0.5d0*abs(d-sqrt(abs(d^2-4*a2*(bN)^2))));

beta1 = sign(b1s*n1+b2s*n2+b3s*n3);

if ( abs(cf*cf-cs*cs) < 1e-12)
    alphaf = sin(atan(1));
    alphas = cos(atan(1));
else
    alphaf = sqrt(abs(a2 - cs*cs))/sqrt(abs(cf*cf-cs*cs));
    alphas = sqrt(abs(cf*cf - a2))/sqrt(abs(cf*cf-cs*cs));
end

TxN1=n3*t2-n2*t3;
TxN2=n1*t3-n3*t1;
TxN3=n2*t1-n1*t2;
%
%    # 1 - right eigenvector Entropy Wave
ru(1,1)   =sqrt(g1/gamma)*rhosq;
ru(2,1)   =0;
ru(3,1)   =0;
ru(4,1)   =0;
ru(5,1)   =0;
ru(6,1)   =0;
ru(7,1)   =0;
ru(8,1)   =0;
%
%    # 2 - right eigenvector Divergence Wave
ru(1,2)   =0;
ru(2,2)   =0;
ru(3,2)   =0;
ru(4,2)   =0;
ru(5,2)   =0;
ru(6,2)   =sq1og*a*n1;
ru(7,2)   =sq1og*a*n2;
ru(8,2)   =sq1og*a*n3;
%
%    # 3 - right eigenvector Alfven Wave
%    lambda = V \times n + b \times n
ru(1,3)   =0;
ru(2,3)   =-sq12*(sqpr*(TxN1));
ru(3,3)   =-sq12*(sqpr*(TxN2));
ru(4,3)   =-sq12*(sqpr*(TxN3));
ru(5,3)   =0;
ru(6,3)   =sq12*sqpor*(TxN1);
ru(7,3)   =sq12*sqpor*(TxN2);
ru(8,3)   =sq12*sqpor*(TxN3);

%
%    # 4 - right eigenvector Alfven Wave
%    lambda = V \times n - b \times n
ru(1,4)   =ru(1,3);
ru(2,4)   =-ru(2,3);
ru(3,4)   =-ru(3,3);
ru(4,4)   =-ru(4,3);
ru(5,4)   =ru(5,3);
ru(6,4)   =ru(6,3);
ru(7,4)   =ru(7,3);
ru(8,4)   =ru(8,3);

%
%    # 5 - right eigenvector
%    lambda = V \times n + C_f
bst = (b1s*t1+b2s*t2+b3s*t3);

ru(1,5)   =sqg2*alphaf*rhosq;
ru(2,5)   =sqg2*((alphaf*a2*n1+alphas*a*((bst)*n1-(bN)*t1)))/(rhosq*cf);

ru(3,5)   =sqg2*((alphaf*a2*n2+alphas*a*((bst)*n2-(bN)*t2)))/(rhosq*cf);

ru(4,5)   =sqg2*((alphaf*a2*n3+alphas*a*((bst)*n3-(bN)*t3)))/(rhosq*cf);
ru(5,5)   =sqg2*alphaf*rhosq*a2;
ru(6,5)   =sqg2*alphas*a*t1;
ru(7,5)   =sqg2*alphas*a*t2;
ru(8,5)   =sqg2*alphas*a*t3;

%
%    # 6 - right eigenvector
%    lambda = V \times n - C_f
ru(1,6) = ru(1,5);
ru(2,6) = -ru(2,5);
ru(3,6) = -ru(3,5);
ru(4,6) = -ru(4,5);
ru(5,6) = ru(5,5);
ru(6,6) = ru(6,5);
ru(7,6) = ru(7,5);
ru(8,6) = ru(8,5);

%
%    # 7 - right eigenvector
%    lambda = V \times n + C_s
ru(1,7) =sqg2*alphas*rhosq;
ru(2,7) =beta1*sqg2*(alphaf*cf^2*t1+a*n1*alphas*(bN))/(rhosq*cf);
ru(3,7) =beta1*sqg2*(alphaf*cf^2*t2+alphas*a*(bN)*n2)/(rhosq*cf);
ru(4,7) =beta1*sqg2*(alphaf*cf^2*t3+alphas*a*(bN)*n3)/(rhosq*cf);
ru(5,7) =a^2*sqg2*alphas*rhosq;
ru(6,7) =-sqg2*alphaf*a*t1;
ru(7,7) =-sqg2*alphaf*a*t2;
ru(8,7) =-sqg2*alphaf*a*t3;
%
%    # 8 - right eigenvector
%    lambda = V \times n - C_s
ru(1,8) =ru(1,7);
ru(2,8) =-ru(2,7);
ru(3,8) =-ru(3,7);
ru(4,8) =-ru(4,7);
ru(5,8) =ru(5,7);
ru(6,8) =ru(6,7);
ru(7,8) =ru(7,7);
ru(8,8) =ru(8,7);

%-   ------------------------------------------------------------------------
%    CONSERVATIVE E_VECTORS
%-   ------------------------------------------------------------------------
%
for m=1:8
    
    R(1,m)=ru(1,m)/g1;
    R(2,m)=(ru(1,m)*u1 + ru(2,m)*rho)/g1;
    R(3,m)=(ru(1,m)*u2 + ru(3,m)*rho)/g1;
    R(4,m)=(ru(1,m)*u3 + ru(4,m)*rho)/g1;
    R(5,m)=(ru(5,m)/g1+B1*ru(6,m)+B2*ru(7,m)+B3*ru(8,m)+0.5d0*ru(1,m)*(u1^2+u2^2+u3^2) + ru(2,m)*u1*rho + ru(3,m)*u2*rho+ru(4,m)*u3*rho)/g1;
    R(6,m)=ru(6,m)/g1;
    R(7,m)=ru(7,m)/g1;
    R(8,m)=ru(8,m)/g1;
    
end

end

function [t1,t2,t3] = compute_t_vector(n1, n2, n3, B1, B2, B3)

    % compute a vector t such that |t|=1, t.n = 0, and t is in span{n, B}

    [nB1,nB2,nB3] = compute_cross(n1, n2, n3, B1, B2, B3);

    if (nB1^2 + nB2^2 + nB3^2 < eps)
       [tp1,tp2,tp3] = pick_orthogonal(n1, n2, n3);
    else
       [tp1,tp2,tp3] = compute_cross(nB1, nB2, nB3, n1, n2, n3);
    end

    tpnorm = sqrt(tp1^2 + tp2^2 + tp3^2);
    t1 = tp1/tpnorm;
    t2 = tp2/tpnorm;
    t3 = tp3/tpnorm;
end
    
function [c1,c2,c3] = compute_cross(a1, a2, a3, b1, b2, b3)
    % compute cross product of a and b, and store result in c

    c1 = a2*b3 - a3*b2;
    c2 = a3*b1 - a1*b3;
    c3 = a1*b2 - a2*b1;
end

        function [r1,r2,r3] = pick_orthogonal(v1, v2, v3)
    % return some r such that v.r = 0
    % assumption:  norm of v is not very small
    if(abs(v1) >= abs(v2)&&abs(v1) >= abs(v3))
       r2 = 1.0;
       r3 = 1.0;
       r1 = (-v2-v3) / v1;
    elseif(abs(v2) >= abs(v1)&&abs(v2) >= abs(v3))
       r1 = 1.0;
       r3 = 1.0;
       r2 = (-v1-v3) / v2;
    elseif(abs(v3) >= abs(v1)&&abs(v3) >= abs(v1))
       r1 = 1.0;
       r2 = 1.0;
       r3 = (-v1-v2) / v3;
    end
    end