% Fisher Memory Matrix and Eigenspectrum
% Linear Randomnly Connected Homogenous Network with Static Synapses 

%%************************************************************************

function [J_x] = DynamicUniform(re_o, Rates, UE, XE)

% Parameters:
tau_d = 0.2; %ms
tau_f = 0.15; %ms
N = 100; % network size;
U = zeros(1,N) + 0.20;
U = diag(U);   
I = eye(N);

%External Input
S_t = zeros(N,1);
S_t(5) = 2;
So = 10; %mv (mean current)
Input = zeros(N,N) + So;
Input(:,1) = Input(:,1) + S_t;
Sigma = zeros(N,N) + randn(N,N);
V = eye(3*N); %Input connectivity 


%Connectivity
meanw = 0; 
variancew = 4;
d = 0.10;
W  = sprandn (N,N,d)*(variancew^1/2) + meanw;
We = zeros(N,N) + W/N;


      
    %Steady States
   
    ue_o = U*(1+tau_f*re_o/1+U*re_o*tau_f);
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
    
    a1 = -I + We*Ds_o;
  
    a2 = We*(Dd_o);

    a3 = We*(Df_o); 
    
    b1 = U.*ue_o;

    b2 = -1/tau_f-U.*re_o;

    b3 = zeros(N,N);

    c1 =Ds_o;
  
    c2 = Dd_o ;

    c3 = -1/tau_d+U*(Ds_o);
    
    %Effective Connectivity Matrices

    J_x = [ a1 a2 a3; b1 b2 b3; c1 c2 c3];
    
    
J_I = zeros(3*N,3*N);
J_I(1,:) = diag(V);

%Effective Time-varying Variables
L = size(Input,2);
deltax = Rates(:,1) - re_o; %deviations at the time of perturbation
deltaUE = UE(:,1) - ue_o;
deltaXE = XE(:,1) - xe_o;
I_o = zeros(N,N) + So;


%Trajectories of the Linearized System
tau_m = 0.006;
dt = 0.0001;
alpha = dt/tau_m;
DeltaX(:,1) = vertcat(deltax, deltaUE, deltaXE);
t(1) = 0;

for n = 1:L

    SynapticInput = zeros(N,1);
    
   
DeltaI(:,n) = vertcat (I_o(:,n), SynapticInput, SynapticInput);
%Input(:,n) - Sigma(:,n) - I_o(:,n);

    t(n+1) = t(n) + dt;
    DeltaX(:,n+1) = DeltaX(:,n) + alpha*(J_x*DeltaX(:,n) + J_I*DeltaI(:,n)*sqrt(dt));  
    DeltaX(DeltaX < 0) = 0; 
    
     % save network activition
     Linear_Re= DeltaX; 
   
end

figure(1)
plot(t,Linear_Re)

%Fisher Information

Cov = zeros(N,N);
    
for k = 1:L
      J_x_k = (J_x)^k;
      
       J_x_K (:,:,k) = J_x_k;
      
      Cov = J_x_K(:,:,k)*(J_x_K(:,:,k))';
   
      Cov_matrix = Cov + Cov(k);
      
end

Covn = inv(Cov_matrix);

for p = 1:L
    
    for l = 1:L
        
        FMM = J_I'*J_x_K(:,:,p)'*(Covn)*J_x_K(:,:,l) *J_I;

      FMC = diag(FMM);
      figure(2)
      plot(FMC)
    end
    
end

     
    
    %plot evalues of J 
     
    evalues = eig(J_x);    % Get the eigenvalues of J

   figure(3)
   plot(real(evalues),imag(evalues),'r*') %   Plot real and imaginary parts
   xlabel('Real')
   ylabel('Imaginary')
     
   figure(4)
   plot(evalues/sqrt(N)*variancew^1/2,'r*') 
   axis([-1.1 1.1 -1.1 1.1])
     
 
end
