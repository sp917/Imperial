function T = profile4(T0,be,lam,z0,h)

%z = -[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1800 2000 2200 2400 2600 2700 2800 3000 3200 3400 3600 3800 4000 4200 4400 4600];
z = -[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600];

zmid = zeros(length(z)-1,1);

for i=1:length(zmid)
    zmid(i) = (z(i) + z(i+1))/2;
end

T00 = 22;
be0 = 19;
lam0 = 300;
h0 = 2500;
z00 = 0;

H = -z(length(z));

T = Tz2(zmid,T0,be,lam,z0) + zmid/h;

T1 = Tz2(zmid,T00,be0,lam0,z00) + zmid/h0;

hold off
plot(T,zmid)
hold on
plot(T1,zmid)




T0,be,z0,lam,h