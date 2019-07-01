z = -[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600]; 

T0 = 22;
al = 19;
L = 300;
h = 2500;
H = 1600;

hinv = 1/h;

T = Tz(z,T0,al,L,h,H);

plot(T,z)
hold on

l = 3;

Lb = 100*l;
hb = 400/l;

hbinv = 1/hb;

alb =  ((H*H/2)*(hbinv-hinv) + al*(L*log(cosh(H/L)) - H*tanh(H/L)))/(Lb*log(cosh(H/Lb)) - H*tanh(H/Lb)) ;

T0b = T0 + alb*tanh(H/Lb) - al*tanh(H/L) + H*(hbinv - hinv);

%{
T0b = T0 + 5;
hb = h/9;
alb = 1.5*al;
Lb = 1.5*L/9;
%}
Tb = Tz(z,T0b,alb,Lb,hb,H);



plot(Tb,z)
hold off

fun = @(z) Tz(z,T0,al,L,h,H);
funb = @(z) Tz(z,T0b,alb,Lb,hb,H);

intT = integral(fun,-H,0);
intTb = integral(funb,-H,0);

intT
intTb

intT = intTz(T0,al,L,h,H)
intTb = intTz(T0b,alb,Lb,hb,H)

T0b
alb
Lb
hb
