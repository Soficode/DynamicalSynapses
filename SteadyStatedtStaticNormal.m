function [re_o, Rates] = SteadyStatedtStaticNormal

% Parameters:
N = 100; % network size;
   
%Initializations
Re = 45*rand(N,1); %Initial Rates

%External Input
S_t = zeros(N,1);
S_t(5) = 200;
So = 10; %mv (mean current)
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t;
Sigma = zeros(N,N)+ randn(N,N);
V = eye(N); %Input connectivity 

%Connectivity

meanw = 0; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = W - tril(W,-1) + tril(W,1)';

%Training
dt =0.0001;
trainTime = size(Input,1);
L = trainTime; 
Rates = zeros(N,L);  % complete state history
t(1) = 0;
Re(:,1) = Re;

for n = 1:L
 
   tau_m = 0.006;
   alpha = 1/tau_m; %ms tau_m
 
   
 % Network activation
 
  h(:,n) = We*Re(:,n) + V*(Input (:,n)) + Sigma(:,n);
 %h(:,n) = sigmf(h(:,n),[2 4]);
  
    t(n+1) = t(n) + dt;
    Re(:,n+1) = Re(:,n) - dt*(alpha*(Re(:,n)-h(:,n)*sqrt(dt)));  
    Re(Re < 0) = 0; 
    
 
     % save network activition
     Rates= Re; 
   
end
 
re_o = Rates(:,L);

plot(t,Re)

   evalues = eig(We);    % Get the eigenvalues of effective connectivity matrix

   figure(3) %   Plot real and imaginary parts
     plot(real(evalues),imag(evalues),'r*') 
     xlabel('Real')
     ylabel('Imaginary')
     
     figure(4) % Check if eigenspectrum lies whithin the unit circle when normalized
     plot(evalues/(sqrt(N)*variancew^1/2),'r*') 
     axis([-1.1 1.1 -1.1 1.1])
   

   
end