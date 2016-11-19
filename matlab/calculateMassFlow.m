function [ p0c, mDot] = calculateMassFlow( rhoOrifice, p0Bottle, T0Combust, AThroat, cp, p0cGuess, cd, AOrifice,F)
%We try and solve for the mass flow and pressure drop in the chamber
%iteratively

%F is the dimensionless mass flow rate @ M=1

err = 1;
tol = 1e-8;
maxIter = 100;

for i = 1:maxIter
    mDot = cd*AOrifice * sqrt(rhoOrifice*2*(p0Bottle - p0cGuess));
    p0c = mDot*sqrt(cp*T0Combust)/(AThroat*F);

    err = p0c - p0cGuess;
    p0cGuess = p0c;
    
    if abs(err)< tol
        break
    end
end

if abs(err) >= tol
    %We ran out of iterations and din't converge
    disp('did not converge')
    error('Convergence issue -bottle pressure too low?')
    
end

