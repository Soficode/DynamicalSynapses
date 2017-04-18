
function [re_o, Rates] = SteadyStateDynamic

% Parameters:
N = 100; % network size;
tau_d = 200; %Synaptic depression time constant
tau_f = 150; %Synaptic facilitation time constant
U = zeros(1,N) + 0.20; %Maximum faciliation parameter
U = diag(U);
   
%Initializations
Re = 45*rand(N,1); %Initial Rates
ue = zeros(N,1); %Initial synaptic depression variables 
xe = zeros(N,1) + 1; %Initial synaptic faciliation variables

%External Input in HZ as in previous static code
S_t = zeros(N,1);
S_t(5) = 2;
So = 10; 
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t;
Sigma = zeros(N,N) + randn(N,N);
V = eye(N);

%Connectivity
meanw = -2; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = zeros(N,N) + W/N;

%Training
dt = 1; % 0.001; 
trainTime = size(Input,1); 
Rates = zeros(N,L);  % complete state history
size(trainTime)
t(1) = 0;
Re(:,1) = Re;
ue(:,1) = ue;
xe(:,1) = xe;

for n = 1:trainTime
 
   tau_m = 60;
   alpha = 1/tau_m; 
    alpha_ue = 1/tau_f; 
     alpha_xe = 1/tau_d;


 % Synaptic variables multiply the weight matrix to form the effective
 % connectivity matrix
   Ds = (ue(:,n).*xe(:,n)); 
   Ds = diag(Ds);
   
 % Network activation
  
 h(:,n) = We*Ds*Re(:,n) + V*Input (:,n)/dt + Sigma(:,n);
 h(:,n) = sigmf(h(:,n),[2 4]);
 
 %Rate Dynamics

    t(n+1) = t(n) + dt;
    Re(:,n+1) = Re(:,n) - dt*(alpha*(Re(:,n)+h(:,n)));  
    Re(Re < 0) = 0; 
    
    
  %Synaptic Dynamics - Short term depression  
    ue(:,n+1) = ue(:,n) + dt*alpha_ue*((- ue(:,n)/tau_f) + U*(1 - ue(:,n)).*Re(:,n));
    ue(ue < 0) = 0;
    
   %Short term facilitation 
    xe(:,n+1) = xe(:,n) + dt*alpha_xe*(((1 - xe(:,n))/tau_d) - xe(:,n).*ue(:,n).*Re(:,n));
    xe(xe > 1) = 1;
    xe(xe < 0) = 0;
 
     % save network activition
     Rates = Re; 
   
end

   
re_o = Rates(:,L);

plot(t,Re)
 
 
end
