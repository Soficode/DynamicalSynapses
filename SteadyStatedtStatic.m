
function [re_o, Rates] = SteadyStatedtStatic

% Parameters:
N = 100; % network size;
   
%Initializations
Re = 45*rand(N,1); %Initial Rates

%External Input in Hz
S_t = zeros(N,1);
S_t(5) = 2; % one current injection at the beginning
So = 10; % Constant background Input
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t; %Input current includes background and transient
Sigma = zeros(N,N)+ randn(N,N) + 200; % External noise
V = eye(N); %Input connectivity 

%Connectivity (Sparse and Random. Weights follow a normal distribution and
%are scaled 1/n)
meanw = 0; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = zeros(N,N) + W/N;

%Training
dt = 0.001; % Integration step size
trainTime = size(Input,1); %length of integration = lenght of the input
L = trainTime; 
Rates = zeros(N,L);  % Initialize matrix that stores history of the rates
t(1) = 0; %time variable
Re(:,1) = Re; %Rates

for n = 1:L
 
   tau_m = 0.060; %membrane time constant
   alpha = 1/tau_m; %leak parameter
 
   
 % Network activation
 
  h(:,n) = We*Re(:,n) + V*(Input (:,n)) + Sigma(:,n); %Input one neuron(External input and recurrent connectivity)
% h(:,n) = sigmf(h(:,n),[2 4]); %Non linear transfer function 
  
    t(n+1) = t(n) + dt; %Time variable with step size dt
    Re(:,n+1) = Re(:,n) - dt*(alpha*(Re(:,n)-h(:,n)*sqrt(dt)));  %Discrete upate of the rates
    Re(Re < 0) = 0; %Rates are always positive
    
 
     % save network activition
     Rates= Re; 
   
end
 
re_o = Rates(:,L); %Save the value of the rates at the last time point - approaches steadystate

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