clc
clear all

start_freq = 6;
end_freq = 35;
Solution = zeros(10001,end_freq - start_freq+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       Swipe frequency
sigma = 0.1;

weighted = [0.002218196, 0.008773135, 0.027023158, 0.064825185,...
    0.12110939, 0.176213123, 0.199675627, 0.176213123, 0.12110939,...
    0.064825185, 0.027023158, 0.008773135, 0.002218196];


for i = 1:13
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       PARAMETERS
    alpha1 = 0.0108;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    1
    alpha2 = 0.0053;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    2
    
    gamma1 = 1.76 * 1e2; %GHz%%%%%%%%%%%%%%%%%%%%%%%    3
    gamma2 = 1.76 * 1e2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%    4
    
    mu0 = 4 * pi * 1e-7;%%%%%%%%%%%%%%%%%%%%%%%%%%%%    5
    
    J = -0.27 * 1e-3;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    6
    
    Hk10 = 0.779 / mu0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    7
    Ms1 = 0.6 / mu0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    8
    Hk1 = Hk10 - Ms1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    9
    d1 = 5.14 * 1e-9;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    10
    
    Hk20mean = -0.12843 / mu0;%%%%%%%%%%%%%%%%%%%%%%    11
    Hk20 = Hk20mean + Hk20mean * sigma * (0.5*i-3.5);
    Ms2 = 1.69 / mu0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    12
    Hk2 = Hk20 - Ms2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    13
    d2 = 3 * 1e-9;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    14
    
    eta1 = J ./ (mu0 * Ms1 * d1);%%%%%%%%%%%%%%%%%%%    15
    eta2 = J ./ (mu0 * Ms2 * d2);%%%%%%%%%%%%%%%%%%%    16
    Hk20 = -0.16843 / mu0;
    parameters = [alpha1, alpha2, gamma1, gamma2, mu0, J, Hk10, Ms1, Hk1, d1,...
        Hk20, Ms2, Hk2, d2, eta1, eta2];
    
    [beta1, beta2] = FindBeta(parameters);
    temp = FindIntensity(parameters, beta1, beta2, start_freq,end_freq);
    Solution = Solution + temp .* weighted(i);
    toc
end

function [beta1, beta2] = FindBeta(parameters)
mu0 = parameters(5);
J = parameters(6);
Ms1 = parameters(8);
Hk1 = parameters(9);
d1 = parameters(10);
Ms2 = parameters(12);
Hk2 = parameters(13);
d2 = parameters(14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       DIFINE FUNCTION
%%% x(1)   x(2)   x(3)   x(4)
%%% phi1   phi2   beta1  beta2
f = @(x,h) -Ms1.*d1.*h.*cos(x(3)).*cos(x(1)) - 0.5.*mu0.*Ms1.*d1.*Hk1.*sin(x(3))^2 ...
    -J.*(cos(x(3)) .* cos(x(1)) .* cos(x(4)) .* cos(x(2))...
    + cos(x(3)) .* sin(x(1)) .* cos(x(4)) .* sin(x(2)) + sin(x(3)) .* sin(x(4)))...
    -Ms2.*d2.*h.*cos(x(4)).*cos(x(2)) - 0.5.*mu0.*Ms2.*d2.*Hk2.*sin(x(4))^2;
options = optimset('MaxFunEvals',10000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       COMPUTE
%%%%%%%%%%%%%%%%%%%%%%%     Angle 1

Numberofpoints = 10001;                                         % MAY CHANGE
field = linspace(0,1,Numberofpoints);
x0 = [0, 0, 0, 0]; % initial finding location
solution = [];
for h = field
    fun = @(x)f(x,h);
    x = fminsearch(fun, x0, options);
    if x(1) > pi / 2
        x(1) = x(1) - pi;
        x(3) = pi - x(3);
        if x(3) > 2 * pi
            x(3) = x(3) - 2 * pi;
        elseif x(3) < -2 * pi
            x(3) = x(3) + 2 * pi;
        end
    end
    if x(2) > pi / 2
        x(2) = x(2) - pi;
        x(4) = pi - x(4);
        if x(4) > 2 * pi
            x(4) = x(4) - 2 * pi;
        elseif x(4) < -2 * pi
            x(4) = x(4) + 2 * pi;
        end
    end
    x0 = x;
    solution = [solution; [x0,fun(x0)]];
end
phi1 = solution(:,1);
phi2 = solution(:,2);
beta1 = solution(:,3);
beta2 = solution(:,4);

end
function Solution = FindIntensity(parameters, beta1, beta2, first, last)
Solution = [];
alpha1 = parameters(1);
alpha2 = parameters(2);
gamma1 = parameters(3);
gamma2 = parameters(4);
mu0 = parameters(5);
Hk1 = parameters(9);
Hk2 = parameters(13);
eta1 = parameters(15);
eta2 = parameters(16);
for k = first:last
    w = 2 * pi * k;
    OMEGA = 1i * w * eye(4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       MATRICS
    GAMMA = [gamma1/(1 + alpha1^2),0,0,0;0,gamma1/(1 + alpha1^2),0,0;...
        0,0,gamma2/(1 + alpha2^2),0;0,0,0,gamma2/(1 + alpha2^2)];
    
    LAMBDA = [1,-alpha1,0,0;alpha1,1,0,0;...
        0,0,1,-alpha2;0,0,alpha2,1;];
    
    stepsize = 1/(size(beta1,1) - 1);
    solution = zeros(size(beta1));
    
    for i = 1 : size(beta1,1)
        he = (i - 1) * stepsize;
        
        c1 = cos(beta1(i));
        c2 = cos(beta2(i));
        c12 = cos(beta1(i) - beta2(i));
        XIe = [0,-he*c1,0,0;he*c1,0,0,0;0,0,0,-he*c2;0,0,he*c2,0;];
        XIa = mu0 * [0,Hk1*cos(2*beta1(i)),0,0;Hk1*sin(beta1(i)).^2,0,0,0;...
            0,0,0,Hk2*cos(2*beta2(i));0,0,Hk2*sin(beta2(i)).^2,0;];
        XIc = mu0 * [0,-eta1*c12,0,eta1*c12;eta1*c12,0,-eta1,0;0,eta2*c12,0,-eta2*c12;...
            -eta2,0,eta2*c12,0;];
        XI = XIa + XIc + XIe;
        solution(i) = 1./abs(det(GAMMA*LAMBDA*XI-OMEGA));
    end
    solution = solution ./ max(solution);%%%%%%%%%%%%%%
    Solution = [Solution solution];
end
end

% % alpha1 = parameters(1);
% % alpha2 = parameters(2);
% % gamma1 = parameters(3);
% % gamma2 = parameters(4);
% % mu0 = parameters(5);
% % J = parameters(6);
% % Hk10 = parameters(7);
% % Ms1 = parameters(8);
% % Hk1 = parameters(9);
% % d1 = parameters(10);
% % Hk20 = parameters(11);
% % Ms2 = parameters(12);
% % Hk2 = parameters(13);
% % d2 = parameters(14);
% % eta1 = parameters(15);
% % eta2 = parameters(16);