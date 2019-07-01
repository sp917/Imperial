function intT = intTz(T0,al,L,h,H)

%The following is the integral of Tz from -H to 0:

intT = T0*H - L*al*log(cosh(H/L)) - H*H/(2*h);

end
