% Fisher Memory Matrix and Eigenspectrum
% Linear Randomnly Connected Homogenous Network with Static Synapses 

%%************************************************************************

function [We] = StaticRandom (Rates,re_o)

% Parameters:
N = 100; 
I = eye(N);

%Connectivity
meanw = 0; 
variancew = 4;
d = 0.10;  
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = zeros(N,N) + W/N;

%External Input
S_t = zeros(N,1);
S_t(5) = 2;
So = 10; 
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t;
Sigma = zeros(N,N) + randn(N,N);
V = eye(N); 

%Linearized System. Rate equation was linearized around the fixed point Ro.
    
%Effective Connectivity Matrices - jacobian matrices with respect to
%perturbations and evaluated around the fixed point

J_x = -I + We*I; %derivative with respect to deviations from Ro (effective connectivity matrix)
J_I = zeros(N,N);
J_I(1,:) = diag(V); %derivative with respect to input perturbation (effective input connectivity)

%Effective Time-varying Variables
L = size(Input,2); %Length of Integration
deltax = Rates(:,1) - re_o; %deviations from Ro at the time of perturbation
I_o = zeros(N,N) + So; %Steady state of the Input


%Trajectories of the Linearized System - evolution of the deviations of the
%rate 
dt = 1; %0.0001;
DeltaX(:,1) = deltax; %Dynamic variable that captures deviations from steady state given a perturbation 
t(1) = 0;

for n = 1:L
 
    
DeltaI(:,n) = I_o(:,n)-Input(:,n)- Sigma(:,n); %Input deviations from steady state - deterministic
size(DeltaX(:,n))
size(J_x)
size(J_I)
size(DeltaI(:,n))

%Dynamics of the linearized system using jacobians as effective weight
%connectivity and input connectivity matrices
    t(n+1) = t(n) + dt;
    DeltaX(:,n+1) = DeltaX(:,n) + dt*(J_x*DeltaX(:,n) + J_I*DeltaI(:,n));  
    DeltaX(DeltaX < 0) = 0; 
    
     % save network activition
     Linear_Re= DeltaX; 
   
end

figure(1)
plot(t,Linear_Re)

%Fisher Information - depends on the connectivity matrix, covariance matrix
%and the input connectivity

Cov = zeros(N,N); %Initialize covariance matrix : Cov = sum(J^k*J^k'). The sum goes from k = 0 to k = infinity - here we use the end of the input(L)
    
for k = 1:L
      J_x_k = (J_x)^k; % Calculates J^k for every k 
      
      Cov = J_x_k*(J_x_k)';
      
      Cov_matrix = Cov + Cov(k);
      
      J_x_K (:,:,k) = J_x_k; %Stores all values of J_x^k
      
end

Covn = inv(Cov_matrix);

%Fisher Information Matrix, where J_I'*J_x^k*Covn*J_x^l*J_I is only entry
%of the matrix at (k,l). To compute the matrix, use the above formula for
%all time steps (all values of k and l for 1:L)

for p = 1:L
    
    for l = 1:L
        
        FMM = J_I'*J_x_K(:,:,p)'*(Covn)*J_x_K(:,:,l) *J_I; %Fisher Memory Matrix
    end
    
 end
    
    FMC = diag(FMM); %The diagonal elements of the matrix are the Fisher Memory Curve
      figure(2)
      plot(FMC)
    

 %plot evalues of J 
     
    evalues = eig(J_x);    % Get the eigenvalues of effective connectivity matrix

   figure(3) %   Plot real and imaginary parts
     plot(real(evalues),imag(evalues),'r*') 
     xlabel('Real')
     ylabel('Imaginary')
     
     figure(4) % Check if eigenspectrum lies whithin the unit circle when normalized
     plot(evalues/(sqrt(N)*variancew^1/2),'r*') 
     axis([-1.1 1.1 -1.1 1.1])


end
