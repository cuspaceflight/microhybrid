%MICRO HYBRID
%V 0.1
%
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% There are a lot of assumptions in this version (0.1). These range from 
% being okay to being pretty terrible
%
%
% Some of the worst assumptions are:
%
% Properties of nitrous: currently we are approximating nitrous with the
% thermodynamic properties of air. These should be updated with correct
% values
%
% Nitrous combustion: Preheat temperature should be set to that of the
% temperature at which nitrous decomposes. Correct energy release should be
% calculated 
%
%
%
%Bottle properties: temperatue of the bottle (nitrous cannister) should
%change during the burn
%
%Nozzle: when calcualting the force there may be shocks in the nozzle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%Nitrous properties%%%%%%%%%%
% TODO - these are properties of air....
cpNitrous = 995; %From http://catalog.conveyorspneumatic.com/
RNitrous = 287; %J/kgK
TNitrous = 300; %K

%%%%%%%%Combustion properties%%%%%%%
T0Preheat = 600; %K
QRelease = 100; %J/kg, energy released when nitrous decomposes
T0Combust = QRelease/cpNitrous + T0Preheat;
cpCombust = 1500; %Typ val of cp for combustion gases (TODO: add chemistry)
gammaCombust = 1.333; %Typ val of gamma for combustion gases

%%%%%Injector discharge coefficient%%%
cd = 0.65; %Taken from Sutton's Rocket Propulsion elements for sharp edged orifice w/ d < 2.5mm
dOrifice = 0.0005;
AOrifice = pi*(dOrifice^2)*0.25; %m^2

%%%%Nozzle properties%%%%%%
rThroat = 0.002; %m
rExit = 0.003; %m
AThroat = pi*rThroat^2; %m^2
AExit   = pi*rExit^2; %m^2
AeOnAt = AExit/AThroat;

%%%%%%%%%%Starting properties%%%%%%%
mNitrous = 0.016; %kg
pNitrous = 6e5; %Pa
bottleVolume = 1e-5; %m^3


%%%%%%%%%%Time stepping properties%%%%%
deltaT = 0.0001;
massStopThreshold = 1e-3*mNitrous;

pBottleHist = [];
mHist = [];
mDotHist = [];
tHist = [];
pChamberHist = [];
FMomHist =[];
FPreHist = [];
FHist = [];
preHeatPower = []; %Power required for preheating

t =0.;

%sonic p threshold
pAmbient = 1e5;
pThreshold = 0.5398; %TODO: Change so it automatically calculates it for given gamma


%Starting p0c guess:
p0c = (pAmbient/pThreshold) * 2.;

%We precalcualte exit properties...
propertiesExit = compressible(AeOnAt,1,'AA',gammaCombust); 
%We calculate the required nondim mass flow at throat for combustion gases
nonDimMassFT = ((gammaCombust)/sqrt(gammaCombust-1))*(1 + 0.5*(gammaCombust-1))^(-0.5*(gammaCombust+1)/(gammaCombust-1));


while (mNitrous > massStopThreshold) && (pAmbient/p0c < pThreshold)
    mDot = 0.;
    
    rhoOrifice = mNitrous/bottleVolume;
    
    pNitrous = rhoOrifice*RNitrous*TNitrous;
    
    pBottleHist = [pBottleHist pNitrous];
    pChamberHist = [pChamberHist p0c];
    mHist = [mHist mNitrous];
    mDotHist = [mDotHist mDot];
    tHist = [ tHist t];
    
    
    [p0c, mDot] = calculateMassFlow(rhoOrifice, pNitrous, T0Combust, AThroat, cpCombust,p0c, cd, AOrifice,nonDimMassFT); %Todo- change 1.281 to the correct value for gamma of combustion gases
    
    %We caluclate the thrust
    
    
 
    pExit = p0c/propertiesExit(2);
    TExit = T0Combust/propertiesExit(4);
    rhoExit = pExit/(TExit*RNitrous);
    %momentum contribution
    vExit = mDot/(AExit*rhoExit);
    FMomentum = vExit*mDot;
    FPressure = (pExit-pAmbient)*AExit;
    Force = FMomentum+FPressure;
    FMomHist = [FMomHist FMomentum];
    FPreHist = [FPreHist FPressure];
    FHist = [FHist Force];
    
    
    %We calculate the instantaneous preheat power needed
    heatingPower = mDot*cpNitrous*(T0Preheat- TNitrous);
    preHeatPower = [preHeatPower heatingPower];
    %We update
    deltaM = mDot*deltaT;
    mNitrous = mNitrous - deltaM;
    
    t = t + deltaT;
    
end

plot(tHist,FMomHist,'r',tHist,FPreHist,'k',tHist,FHist,'b')
legend('Momentum contribution','Pressure contribution','Total')
xlabel('Time (s)')
ylabel('Force (N)')
title('Thrust')

figure
plot(tHist,preHeatPower,'k')
xlabel('Time (s)')
ylabel('Power (W)')
title('Instantaneous preheat power requirement (perfect heat transfer)')

%We caluclate the impulse
impulse = trapz(tHist,FHist);

disp('Umpulse is (NS):');
disp(impulse);