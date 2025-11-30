function [obj] = Objective_BM(t,x,u,par)
% [obj, objGrad] = pathObjective(u)
%
% Computes the objective function (and gradients) for the simple pendulum
%
len = length(x(1,:));

c = linspace(0,1,len);

obj = -(VV(x(2,:), par) - ( x(1,:) + t .* x(2,:) ) + par.alpha*x(1,:) - par.K * (x(2,:)>0) ) .* ff(c,par);


end
