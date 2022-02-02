%% Parenting as investment (Cunha 2007), Yifan Lyu, SSE
clear; clc;
warning off;

% define parameters
lambda = 0.5; % risk aversion in utility function
w = 1;          % wage rate
%theta_p = 2;    % parent's human capital
b = 1;          % initial wealth
b_prime = 1;    % wealth left for child
r = 0.05;       % interest rate
beta = 0.96;    % discount factor
v = 0.5;        % parental altrusim
g = 0.1;        % rate of wage increase

% human capital production function
gamma = 0.3*ones(3,2); % skill multiplier
phi = 0.3;       % substitution parameter, asssume constant

%theta_1 = 1;     % initial draw of human capital
theta_2 = @(I_1, theta_p, theta_1) ...
       (gamma(1,1)*theta_1'.^phi + gamma(2,1)*I_1.^phi + gamma(3,1)*theta_p.^phi).^(1/phi);
theta_3 = @(I_1, I_2, theta_p, theta_1) ...
       (gamma(1,2)*theta_2(I_1, theta_p,theta_1).^phi + gamma(2,2)*I_2.^phi + gamma(3,2)*theta_p.^phi).^(1/phi);

% set up maximisation problem
c1 = @(I_1,a,theta_p) w*theta_p + b - I_1 - a/(1+r);
c2 = @(c1,theta_p) c1.*( 1/((1+r)*beta) ).^(1/(lambda-1)); % from Euler equation
I2 = @(c2,a,theta_p)    w*(1+g)*theta_p + a - c2 - b/(1+r);

% guess a form of value function: increasing in both argument
% create grid for theta_p and theta_1:
N_g = 20;  % number of grid point
theta_p = linspace(0.01,5,N_g); % grid of parental endowment
theta_1 = linspace(1,5,N_g); % grid of new distribution
V_guess = theta_p + theta_1'; % first guess of value function
dis = -inf; % initial distance between V_guess and V_update


for loop = 1
theta_p = linspace(0.01,5,N_g); % grid of parental endowment
theta_1 = linspace(1,5,N_g); % grid of new distribution
I1_guess = linspace(0,1,N_g);
a_guess  = linspace(1,3,N_g);

c1_val = c1(I1_guess,a_guess,theta_p); % guessed value of c1
c2_val = c2(c1_val,theta_p); % guessed value of c2
I2_val = I2(c2_val,a_guess,theta_p); % guessed value of I2
theta_3_val = theta_3(I1_guess,  I2_val,  theta_p, theta_1); 
% guessed value of parental endowment: theta_3_val is a matrix, column is
% theta_1 value
%A = [min(theta_3_val,[],'all'),max(theta_3_val,[],'all')];
% create interpolation
[X,Y] = ndgrid(theta_p, theta_1);
% create griddled interpolant
findVq = griddedInterpolant(X,Y,V_guess);

% create quiry points
%[Xq,Yq] = meshgrid(theta_3_val(:),theta_1);


i = 0;
j = 0;
V_guess = nan(N_g,N_g);
for theta_p = linspace(0.01,5,N_g) % range of parental endowment
    i = i+1;
    for theta_1 = linspace(1,5,N_g) % theta 1 grid given
    
    if j<N_g
    j = j+1;
    else
    j = 1;  
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_guess(j,i) = theta_p + theta_1;
    ini_guess = [I1_guess(i), a_guess(i)]; % initial guess of I1 and a scalar
    
    %c1_val = c1(guess(1),guess(2),theta_p); % guessed value of c1 scalar
    %c2_val = c2(c1_val,theta_p); % guessed value of c2 scalar
    %I2_val = I2(c2_val,guess(2),theta_p); % guessed value of I2 scalar

    %V_update = util(c1_val) + beta*util(c2_val) ...
    %     + beta^2*v*findVq(theta_p,theta_1); % scalar of updated value function
    V_obj = @(guess) - (  util(c1(guess(1),guess(2),theta_p)) + beta*util(c2(c1(guess(1),guess(2),theta_p),theta_p)) ...
             + beta^2*v*findVq(theta_p,theta_1)  );
    options = optimset('display','none');
    [Output,V_output] = fmincon(V_obj,ini_guess, [1,(1+r)], w*theta_p + b,[],[],[0,0],[],[],options);
    V_update(i,j) = -V_output;
    I1_store(i,j) = Output(1);
    a_store(i,j)  = Output(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

% find distance
dis = max( abs(V_guess-V_update) ,[],'all');
fprintf('current distance = %5f \n', dis);
V_guess = V_update;

end



function val = util(c)
logic = (c>0);
lambda = 0.5;
val(logic==1) = (c(logic==1).^lambda-1)/lambda; % parental utility function
val(logic==0) = -1e10;
end