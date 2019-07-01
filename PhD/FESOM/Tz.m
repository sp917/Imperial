function T = Tz(z,T0,al,L,h,H)

T = T0 + al*tanh(z/L) + z/h;

%szL = sech(z/L);
%dT = al*(szL.*szL)/L + 1/h;

%The following is the integral of Tz from -H to 0:

intT = T0*H - L*al*log(cosh(H/L)) - H*H/(2*h);

end

