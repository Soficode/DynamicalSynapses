%% Random Matrix : Linearized Dynamic Synapses Model
%%************************************************************************

function [J_x, DeltaX] = RandomOrthogonal (re_o, Re)

% Parameters:
tau_d = 0.200; %ms
tau_f = 0.15; %ms
tau_m = 0.006; %ms?
N = 1000; % network size;
U = zeros(N,1) + 0.20;   
X = zeros(N,1); 
I = eye(N);

%External Input
S_t = zeros(N,1);
S_t(5) = 20;
S_t(20) = 20;
S_t(50) = 20;
S_t(500) = 20;
S_t(250) = 20;
S_t(340) = 20;
So = 100; % Constant background Input
So = zeros(N,1) + So;
Input = So + S_t; %Input current includes background and transient
epsilon = 0.5;
Sigma = randn(N,1)*epsilon^1/2; % External noise
V = zeros(N,1); %Input connectivity 
V(1) = 1;

%Connectivity
meanw = 0; 
variancew = 4;%0.999/4*N;
d = 0.10;
We  = orth(randn (N,N))*(variancew^1/2) + meanw;


      
    %Steady States
   
    ue_o = U.*(1+tau_f*re_o/1+U.*re_o*tau_f);
    
    xe_o = 1/(1+(ue_o.*re_o*tau_d));
    xe_o = xe_o';
 
    Ds_o = (ue_o.*xe_o);
    Ds_o = diag(Ds_o);
    
    Df_o = (ue_o.*re_o);
    Df_o = diag(Df_o);
   
    Dd_o = (re_o.*xe_o);
    Dd_o = diag(Dd_o);
  
       
     
    %Linearized System - Jacobian Matrix (effective connectivity matrix of
    %the linearized system)
    
    a1 = 1/tau_m*(-I + We*Ds_o);
  
 
    a2 = 1/tau_m*(We*(Dd_o));

    a3 =1/tau_m*( We*(Df_o)); 
    
    b1 = 1/tau_f*(U*ue_o');

    b2 = 1/tau_f*(-1/tau_f-U*re_o');
  
    b3 = 1/tau_f*(zeros(N,N));

    c1 =1/tau_d*(Ds_o);
  
    c2 = 1/tau_d*(Dd_o );

    c3 = 1/tau_d*(-1/tau_d+diag(U)*(Ds_o));
  
    
 %Effective Connectivity Matrices

J_x = [ a1 a2 a3; b1 b2 b3; c1 c2 c3];  
V_s = zeros(N,1);
J_I = vertcat(V, V_s, V_s);

%Effective Time-varying Variables
dt = 0.001; % Integration step size
L = size(Input,1); %length of integration = lenght of the input
deltax = Re(:,1) - re_o; %deviations at the time of perturbation
deltaUE = U - ue_o;
deltaXE = X - xe_o;


%Trajectories of the Linearized System
DeltaX(:,1) = vertcat(deltax, deltaUE, deltaXE);
t(1) = 0;

for n = 1:L
    
Input_delta = So(n) - Input(n) - Sigma(n);   

    t(n+1) = t(n) + dt;
h = J_x*DeltaX(:,n) + J_I*Input_delta;
   
    DeltaX(:,n+1) = DeltaX(:,n) + dt*(h);  
    DeltaX(DeltaX < 0) = 0; 
     
end

figure(1)
plot(t,DeltaX)

    %plot evalues of J 
     
    evalues = eig(J_x);    % Get the eigenvalues of J

   figure(2)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
   ylabel('Imaginary')
     
end