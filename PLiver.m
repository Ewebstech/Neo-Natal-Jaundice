function P = PLiver(t, P0, r_m)
% P0 = Baby's Ability to excrete bilirubin from the liver/bile into the stool at birth is P0
% r_m = Maturity Rate
% t = time
P = P0/(P0 + (1-P0)*exp(-r_m*t));
end
