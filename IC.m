function dx = IC(x,u)
% dx = IC(x,u)
%
% local IC conditions using the envelope theorem
%

pft = x(1,:); % profit, Pi
q = x(2,:); % allocation, q


dq = u;

%dpi = -q


dx = [-q;dq];

end