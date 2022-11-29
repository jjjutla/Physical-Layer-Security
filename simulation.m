close all;
clear all;
R = 5;       %Information rate
x = 2.^R-1;  %Outage thresholds
Omega1 = 1;  %Average channel gain between transmitter 1 and relay
Omega2 = 1;  %Average channel gain between transmitter 2 and relay
N1 = 1;      %Normalized noise variance at transmitter 1
N2 = 1;      %Normalized noise variance at transmitter 2
N3 = 1;      %Normalized noise variance at relay
P1_dB = 0:2.5:50;    %SNR between transmitter 1 and relay
P2_dB = P1_dB;       %SNR between transmitter 2 and relay
P1 = 10.^(P1_dB/10); %Corresponding transmit power range at transmitter 1
P2 = 10.^(P2_dB/10); %Corresponding transmit power range at transmitter 2
P3 = 0.5*P2;         %Relay transmit power (half of transmitters)
K3t = [0 0.05 0.1 0.2]; %Error Vector Magnitudes for transmission at the relay
K3r = K3t; %Error Vector Magnitudes for reception at the relay


%Number of realizations in Monte-Carlo simulations
nbrOfRealizations = 1000000;

%Channel fading realizations
h1 = sqrt(Omega1/2)*(randn(1,nbrOfRealizations) + 1i*randn(1,nbrOfRealizations)); 
h2 = sqrt(Omega2/2)*(randn(1,nbrOfRealizations) + 1i*randn(1,nbrOfRealizations)); 
%Squared norm fading power of h1 and h2
RHO1 = abs(h1).^2;  
RHO2 = abs(h2).^2;  

%Modulation parameters for symbol error rate computation 
alpha = 1; 
beta = 1; 

SNDR1 = zeros(length(P1),nbrOfRealizations,length(K3t)); 
SNDR2 = zeros(length(P1),nbrOfRealizations,length(K3t)); 
PoutSimulation1 = zeros(length(P1),length(K3t)); 
PoutSimulation2 = zeros(length(P1),length(K3t)); 

%Loop all cases of transmit power at the transmitters and assume that
%the power on the relay is varied in the same way. 
for p = 1:length(P1)
    for k = 1:length(K3t)
        
        %Relaying gain as in Eq. (5)
        G = sqrt( P3(p) ./ ( (RHO1*P1(p)+RHO2*P2(p))*(1+K3r(k)^2) + N3 ) );
        
        a1 = (N3/P2(p)) * (1+K3t(k)^2);
        a2 = (N3/P1(p)) * (1+K3t(k)^2);
        b1 = (N1/P3(p)) * (1+K3r(k)^2);
        b2 = (N2/P3(p)) * (1+K3r(k)^2);
        c = K3t(k)^2 + K3r(k)^2 + K3t(k)^2 *K3r(k)^2;
        
        %Compute SNDRs for different fading realizations
        SNDR1(p,:,k) = (RHO1.*RHO2) ./ ( RHO1.^2*(P1(p)/P2(p))*c + RHO1.*RHO2 *c + RHO2*b1 + RHO1*(a1+(P1(p)/P2(p))*b1) + N1*N3/(P2(p)*P3(p)) );
        SNDR2(p,:,k) = (RHO1.*RHO2) ./ ( RHO2.^2*(P2(p)/P1(p))*c + RHO1.*RHO2 *c + RHO1*b2 + RHO2*(a2+(P2(p)/P1(p))*b2) + N2*N3/(P1(p)*P3(p)) );
        
        PoutSimulation1(p,k) = sum(SNDR1(p,:,k) < x) /nbrOfRealizations;
        PoutSimulation2(p,k) = sum(SNDR2(p,:,k) < x) /nbrOfRealizations;    
    end
    
end

figure(1); hold on; box on;

for k = 1:length(K3t)
    plot(P1_dB,PoutSimulation1(:,k),'k*','Markersize',6); %Markers based on Monte Carlo simulations
end

axis([P1_dB(1) P1_dB(end) 0 1.1]);
axis([0 50 -0.02 1.1])
