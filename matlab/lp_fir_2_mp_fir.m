function [hn_min] = lp_fir_2_mp_fir(hn_lin)

% This function helps in converting a linear phase 
% FIR to a minimum phase FIR. It essentially 
% brings all zeros which are outside the unit circle to 
% inside of unit circle.

r = roots(hn_lin);
k = abs(r) > 1.01;
r1 = r(k); % roots outside of unit circle
r2 = r(~k); % roots on or within the unit circle
s = [1./r1; r2];
hn_min = poly(s);
hn_min = hn_min*sqrt(sum(hn_lin.^2))/sqrt(sum(hn_min.^2));