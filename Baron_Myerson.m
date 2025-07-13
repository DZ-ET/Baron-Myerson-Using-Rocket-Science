% Baron_Myerson without fixed cost

%% Initialization

clc; clear;
addpath /Users/dihanzou/Library/CloudStorage/OneDrive-Personal/MATLAB_WD/OptimTraj-master

% market primitives

par.d = 0.001;
par.v = 0:par.d:1;

par.alpha = 1; % welfare weight on profit
% \alpha = 0: strongest redistributional motive; \alpha = 1: utilitarian regulator

%par.CAP = 0.1; % Cap of lump-sum transfer. UNCOMMENT THE PATH CONSTRAINT
%IF NEEDED

% linear demand curve
par.g = ones(1,1001);   % consumer value distribution (pdf)
par.G = par.v;          % consumer value distribution (cdf)


%%% CHOOSE THE COST DISTRIBUTION FROM UNIFORM OR NORMAL

% uniform cost

par.F = par.v;
par.f = ones(1,1001);

%truncated normal

% mu = 0.5; %0.15
% 
% sigma = 0.1;
% 
% xi = (par.v - mu)/sigma;
% beta = (1 - mu)/sigma;
% alpha = -mu/sigma;
% 
%  par.F = (normcdf(xi) - normcdf(alpha))/(normcdf(beta)-normcdf(alpha));
%  par.f = normpdf(xi)/sigma/(normcdf(beta)-normcdf(alpha));


% Johnson-Myatt % not converging
% % 
% par.f = 0.5*normpdf(par.v,0.3,0.1) + 0.5*normpdf(par.v,0.7,0.1);
% 
% par.F = cumtrapz(par.f)*par.d;




%%



% dynamics and objective functions
problem.func.dynamics = @(t,x,u)( IC(x,u) );
problem.func.pathObj = @(t,x,u)( Objective_BM(t,x,u,par) ); % beware of minimization

% problem.func.pathCst = @(t,x,u)( pathConstraint(x,t,par) ); % UNCOMMENT THIS IF THERE IS A PATH CONSTRAINT

% Problem bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 1;
problem.bounds.finalTime.upp = 1;

problem.bounds.state.low = [0; 0];
problem.bounds.state.upp = [inf; 1];


problem.bounds.initialState.low = [0; 0];
problem.bounds.initialState.upp = [20; 1];
problem.bounds.finalState.low = [0;0];
problem.bounds.finalState.upp = [0;0];

problem.bounds.control.low = -inf; % -inf
problem.bounds.control.upp = 0; 

% Guess at the initial trajectory
problem.guess.time = [0,1];
problem.guess.state = [0.2, 0; 0.8, 0];
problem.guess.control = [0, 0];


%%

%%%% Switch between a variety of methods

method = 'trapezoid';
% method = 'trapGrad';   
% method = 'hermiteSimpson';
% method = 'hermiteSimpsonGrad';   
% method = 'chebyshev';   
% method = 'rungeKutta';  
% method = 'rungeKuttaGrad';
% method = 'gpops';

problem.options.defaultAccuracy = 'high';

% % 
% % 
% %% Method-independent options:
problem.options(1).nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-3,...
    'MaxFunEvals',1e4);   %options for fmincon
problem.options(2).nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-6,...
    'MaxFunEvals',5e6);   %options for fmincon




switch method
    
    case 'trapezoid'
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options  % 10
        
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 100;  %method-specific options  % 25
        
    case 'trapGrad'  %trapezoid with analytic gradients
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 45;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 50;  %method-specific options 15
        
    case 'hermiteSimpsonGrad'  %hermite simpson with analytic gradients
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
        
    case 'chebyshev'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'chebyshev'; % Select the transcription method
        problem.options(1).chebyshev.nColPts = 9;  %method-specific options 9
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'chebyshev'; % Select the transcription method
        problem.options(2).chebyshev.nColPts = 30;  %method-specific options 15
        
    case 'multiCheb'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'multiCheb'; % Select the transcription method
        problem.options(1).multiCheb.nColPts = 6;  %method-specific options
        problem.options(1).multiCheb.nSegment = 4;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'multiCheb'; % Select the transcription method
        problem.options(2).multiCheb.nColPts = 9;  %method-specific options
        problem.options(2).multiCheb.nSegment = 4;  %method-specific options
        
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
    
    case 'rungeKuttaGrad'
      
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    case 'gpops'
        problem.options = [];
        problem.options.method = 'gpops';
        problem.options.defaultAccuracy = 'high';
        problem.options.gpops.nlp.solver = 'snopt';  %Set to 'ipopt' if you have GPOPS but not SNOPT
        
    otherwise
        error('Invalid method!');
end




%% Solve the problem
soln = optimTraj(problem);
c = soln(end).grid.time;
pft = soln(end).grid.state(1,:);
q = soln(end).grid.state(2,:);
u = soln(end).grid.control;

cst_val = 1 - GG(pft./q + c, par)  - q;

%% Plot the solution:
figure(1); clf;

subplot(4,1,1)
plot(c,pft)
ylabel('$\Pi$', 'interpreter', 'latex')
title('profit');

subplot(4,1,2)
plot(c,q)
ylabel('$q$', 'interpreter', 'latex')
%ylim([0 0.6])

subplot(4,1,3)
plot(c,u)
ylabel('$\dot{q}$', 'interpreter', 'latex')


% IGNORE THE LAST PANEL IF NO PATH CONSTRAINT
subplot(4,1,4)
plot(c,cst_val)
ylabel('constraint', 'interpreter', 'latex')



