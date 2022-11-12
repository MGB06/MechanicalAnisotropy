clc;clear;close all
 
% load CTAn results
load CTAn_data
 %% Zysset Curnier model with model parameters from Gross et al.
l0s=4609.2;
l1s=3691.6;
m0s=3738;
ks=1.60;
ls=0.99;
 for t=1:length(AreaCheck_CTAn_MIL3)
aval=[m1(t) m2(t) m3(t)];
for w=1:3;
    for j=1:3;
        if w==j
            Stiffness(w,w,t) = (l0s+2*m0s)*rho^ks*aval(w)^(2*ls);
        else
            Stiffness(w,j,t) = l1s*rho^ks*aval(w)^(ls)*aval(j)^(ls);
        end
    end
end
 Stiffness(4,4,t) = 2*m0s*rho^ks*aval(2)^(ls)*aval(3)^(ls);
Stiffness(5,5,t) = 2*m0s*rho^ks*aval(3)^(ls)*aval(1)^(ls);
Stiffness(6,6,t) = 2*m0s*rho^ks*aval(2)^(ls)*aval(2)^(ls);
 clear aval i j
end
