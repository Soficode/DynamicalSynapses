%% Random Symmetric Matrix 

%Choose V as the maximum principal component of Js
%%************************************************************************

function [re_o] = Delay

% Parameters:
N = 1000; % network size;
   
%Initializations
Re = 45*rand(N,1); %Initial Rates

%External Input in Hz
S_t = zeros(N,1);
S_t(5) = 20;
S_t(20) = 20;
S_t(50) = 20;
S_t(500) = 20;
S_t(250) = 20;
S_t(340) = 20;% one current injection at the beginning
So = 100; % Constant background Input
So = zeros(N,1) + So;
Input = So + S_t; %Input current includes background and transient
epsilon = 0.5;
Sigma = randn(N,1)*epsilon^1/2; % External noise
V = zeros(N,1); %Input connectivity 
V(1) = 1;


%Connectivity (Sparse and Random. Weights follow a normal distribution and
%are scaled 1/n)
meanw = 0; 
variancew = 4;%0.999/4*N;
d = 0.10;
We  = orth(randn (N,N))*(variancew^1/2) + meanw;

%Training
dt = 0.001; % Integration step size
L = size(Input,1); %length of integration = lenght of the input
Rates = zeros(N,L);  % Initialize matrix that stores history of the rates
t(1) = 0; %time variable
Re(:,1) = Re;

for n = 1:L
 
   tau_m = 0.060; %membrane time constant
   alpha = 1/tau_m; %leak parameter
 
   
 % Network activation

  h(:,n) = We*Re(:,n) + V*(Input (n)) + Sigma(n); %Input one neuron(External input and recurrent connectivity)
% h(:,n) = sigmf(h(:,n),[2 4]); %Non linear transfer function 
  
    t(n+1) = t(n) + dt; %Time variable with step size dt
   Re(:,n+1) = Re(:,n) - dt*(alpha*(Re(:,n)-h(:,n)*sqrt(dt)));  
    Re(Re < 0) = 0; %Rates are always positive
   
   
end

re_o = Re(:,L); %Save the value of the rates at the last time point - approaches steadystate

figure(1)
plot(t,Re)

   evalues = eig(We);    % Get the eigenvalues of effective connectivity matrix

   figure(2) %   Plot real and imaginary parts
     plot(real(evalues),imag(evalues),'r*') 
     xlabel('Real')
     ylabel('Imaginary')
     
    
end