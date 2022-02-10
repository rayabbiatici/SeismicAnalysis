function [Force,dmat,svars] = ...
   ClassicalPlasticity(Disp,props,fold,svars)

E          = props(1); % Spring Modulus
K          = props(2); % Hardening Modulus
YieldForce = props(3); 
OldPlasticDisp = svars(1);
OldAlpha = svars(2);

TrialForce = E*(Disp - OldPlasticDisp);
ftrial     = abs(TrialForce) - (YieldForce + K*OldAlpha);
      
if ftrial < 0  % elastic step
   Force    = TrialForce; 
   dmat     = E;
   svars(1) = OldPlasticDisp;
   svars(2) = OldAlpha;
   
else           % plastic step
   dgamma   = ftrial / (E+K);
   Force    = (1-dgamma*E/abs(TrialForce))*TrialForce;
   svars(1) = OldPlasticDisp + dgamma*sign(TrialForce);         
   svars(2) = OldAlpha + dgamma;
   dmat     = E*K / (E+K);
   
end


return