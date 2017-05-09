function [re_o, Rates] = SteadyStateDynamicNormal

% Parameters:
N = 100; % network size;
tau_d = 0.200; %ms
tau_f = 0.150; %ms
U = zeros(1,N) + 0.20;
U = diag(U);
   
%Initializations
Re = 45*rand(N,1); %Initial Rates
ue = zeros(N,1);
xe = zeros(N,1) + 1;

%External Input
S_t = zeros(N,1);
S_t(5) = 2;
So = 10; %mv (mean current)
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t;
Sigma = zeros(N,N) + randn(N,N);
V = eye(N); %Input connectivity 

%Connectivity
meanw = 0; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = W - tril(W,-1) + tril(W,1)';


%Training
dt = 0.0001;
trainTime = size(Input,1);
L = trainTime; 
trainTime = size(Input,1);
Rates = zeros(N,L);  % complete state history
size(trainTime)
t(1) = 0;
Re(:,1) = Re;
ue(:,1) = ue;
xe(:,1) = xe;

for n = 1:L
 
   tau_m = 0.006;
   alpha = 1/tau_m; %ms tau_m
    alpha_ue = 1/tau_f; %ms tau_m
     alpha_xe = 1/tau_d; %ms tau_m


   Ds = (ue(:,n).*xe(:,n));
   Ds = diag(Ds);
   
 % Network activation
  
 h(:,n) = We*Ds*Re(:,n) + V*Input (:,n) + Sigma(:,n);
 %h(:,n) = sigmf(h(:,n),[2 4]);

    t(n+1) = t(n) + dt;
    Re(:,n+1) = Re(:,n) - dt*(alpha*(Re(:,n)+h(:,n)*sqrt(dt)));  
    Re(Re < 0) = 0; 
    
    ue(:,n+1) = ue(:,n) + dt*alpha_ue*((- ue(:,n)/tau_f) + U*(1 - ue(:,n)).*Re(:,n));
    ue(ue < 0) = 0;
    
    xe(:,n+1) = xe(:,n) + dt*alpha_xe*(((1 - xe(:,n))/tau_d) - xe(:,n).*ue(:,n).*Re(:,n));
    xe(xe > 1) = 1;
    xe(xe < 0) = 0;
 
     % save network activition
     Rates = Re; 
   
end

   
re_o = Rates(:,L);

plot(t,Re)

   evalues = eig(We*Ds);    % Get the eigenvalues of effective connectivity matrix

   figure(3) %   Plot real and imaginary parts
     plot(real(evalues),imag(evalues),'r*') 
     xlabel('Real')
     ylabel('Imaginary')
     
     figure(4) % Check if eigenspectrum lies whithin the unit circle when normalized
     plot(evalues/(sqrt(N)*variancew^1/2),'r*') 
     axis([-1.1 1.1 -1.1 1.1])
   
 
 
end
