function [c,ceq] = pathConstraint(x,t,par)
% [c] = pathConstraint(x)
% c -- inequality constraint s.t. c <= 0
% ceq -- equality constraint s.t. c == 0

%

pft = x(1,:);
q = x(2,:);

%%% cap of lump-sum transfer

c = pft + t .* q - par.CAP;  

c(q==0) = -1;
ceq = []; 



%%% WEI-ZOU

% c = -1 + GG(pft./q + t, par)  + q;  
% 
% c(q==0) = -1;
% ceq = [];

%%% AMADOR-BAGWELL -- need a good initial guess
% ceq = -1 + GG(pft./q + t, par)  + q;  
% 
% ceq(q==0) = -1;
% 
% c = [];


end