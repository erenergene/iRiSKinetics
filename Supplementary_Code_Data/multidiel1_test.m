addpath(genpath(pwd))

% function [Gamma1,Z1] = multidiel1(n,L,lambda,theta,pol)
clearvars
%%
n = [n_Air(1,2);n_SiO2(1,2);n_Si(1,2)];
L = 0:0.1:300;
lambda = 400:0.1:700;
%%
if nargin==0, 
    help multidiel1;
    return; 
end

if nargin<=4,
    pol='te';
end

if nargin==3,
    theta=0;
end

if size(n,2)==1, 
    n = n';
end 
%%
M = length(n)-2;                                % number of slabs

if M==0,
    L = [];
end                            % single interface, no slabs

%%

theta = theta * pi/180;
                                                                
% costh = conj(sqrt(conj(1 - (n(1) * sin(theta) ./ n).^2)));    % old version 
                                                                
costh = sqrte(1 - (n(1) * sin(theta) ./ n).^2);                 % new version - 9/14/07

if pol=='te' | pol=='TE',
    nT = n .* costh;                            % transverse refractive indices
else
    nT = n ./ costh;                          % TM case, fails at 90 deg for left medium
end

if M>0,
    L = L .* costh(2:M+1);                      % n(i)*l(i)*cos(th(i))
end

r = -diff(nT) ./ (diff(nT) + 2*nT(1:M+1));      % r(i) = (n(i-1)-n(i)) / (n(i-1)+n(i))   

Gamma1 = r(M+1) * ones(1,length(lambda));       % initialize Gamma at right-most interface

for i = M:-1:1,
    delta = 2*pi*L(i)./lambda;                  % phase thickness in i-th layer
    z = exp(-2*j*delta);                          
    Gamma1 = (r(i) + Gamma1.*z) ./ (1 + r(i)*Gamma1.*z);
end

Z1 = (1 + Gamma1) ./ (1 - Gamma1);

% end
