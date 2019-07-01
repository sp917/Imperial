function profile2(z0)

z = -[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600]; 

T0 = 22;
al = 19;
L = 300;
h = 2500;
H = 1600;

hinv = 1/h;

T = Tz(z,T0,al,L,h,H);

hold off
plot(T,z)
hold on

lam = L/9;
be = al*tanh(H/L) + H/h;
be = be/(tanh((z0+H)/lam) - tanh(z0/lam));
%be = 3*al;


T0b = T0 + be*tanh(z0/lam); 

T2 = Tz2(z,T0b,be,lam,z0);

ylim([-1600,0])

T0b
be
lam
z0

plot(T2,z)
hold off

legend('T1','T2')


