% visualise

clear all

%loc = '~/cluster/gold5/';
%loc = '~/';
%loc = '~/Documents/GulfStream/MOM/Code/'
%loc = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/N513/relax/';
loc = '~/OneDrive/Documents/Imperial/PhD/MOM/NoRelax129';

cd(loc);

%outloc = '~/Documents/GulfStream/GOLD/Images/';
%outloc =  '~/folder' %'~/Documents/GulfStream/MOM/Images/';
outloc = '~/OneDrive/Documents/Imperial/PhD/MOM/NoRelax129';

%T = ncread('timestats.nc','Time');
files = dir('prog__*');
nf = size(files);
nf = nf(1);

x = ncread(files(1).name,'xq');
y = ncread(files(1).name,'yq');

%nt = size(T,1);
nx = size(x,1);
ny = size(y,1);
L = 3840000.0; 

% BG state
H1 = 300.0;
f0 = 0.44e-4;
beta = 2e-11;
f = f0 + beta * y;
Q = zeros(ny,nx);
for i = 1:nx
    Q(:,i) = f(:) / H1;
end

q1_scale = Q(1);

%%

loop_start = nf-1;
loop_end = nf-1;
RUN = 1;    % Choose to define PV from a completed RUN (1) or a RUN in progress (2)
if RUN == 1
    PV_av = zeros(nx,ny);
    count = 0;
    %for i = 3:3
    for i = loop_start:loop_end
        disp(i);
        PVnew = ncread(files(i).name,'PV');
        unew = ncread(files(i).name,'u');
        PVnew = PVnew(:,:,1,:);
        unew = unew(:,:,1,:);
        nn = size(PVnew,4);
        for ti = 1:nn
            PV1(:,:,count+ti)=transpose(PVnew(:,:,ti));% - Q(:,:);
            u1(:,:,count+ti)=transpose(unew(:,:,ti));
        end
        count = count + nn;
    end
    nt = size(PV1,3);
end



%%

% Calcualte the PV average, from some time t > 0
PV_av = zeros(ny,nx);
n0 = 1;
for ti = n0:nt
    PV_av = PV_av + PV1(:,:,ti);
end
PV_av = PV_av / (nt - n0 + 1);

%%

ti = nx - floor(nx/10); 
fs = 18;

ts = nt;

% [PV_LS, PV_SS] = SpectralFilter(PV1(1:nx-1,1:ny-1,ts));
% [PV_LS, PV_SS] = BoundaryAdaptFilter(PV1(1:nx-1,1:ny-1,ts),17);
% subplot(131)
% surf(x(1:nx-1),y(1:ny-1),PV_LS,'edgecolor','none'); view(0,90);...
%     shading interp; colorbar(); colormap(jet); axis image;
% subplot(132)
% surf(x(1:nx-1),y(1:ny-1),PV_SS,'edgecolor','none'); view(0,90);...
%     shading interp; colorbar(); colormap(jet); axis image;
% subplot(133)
% surf(x,y,PV1(:,:,ts),'edgecolor','none'); view(0,90);...
%     shading interp; colorbar(); colormap(jet); axis image;
% pause



qlim1 = min(min(PV1(:,:,nt)));
qlim2 = max(max(PV1(:,:,nt)));

figure(1);
surf(x,y,PV1(:,:,ts)/q1_scale,'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y');...
    %caxis([f(1)/H1,5.0e-7]);...
    caxis([1,3.5]);...
    text(x(ti),y(ti),10*PV1(ti,ti,nt)/q1_scale,'$q_{1}$','Interpreter','latex','FontSize',fs,'Color','k');...
    %caxis([qlim1 qlim2]);...
    set(gca,'xTick',x(1):x(nx)-x(2):x(nx-1)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(2):y(ny-1)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'PV_snapshot'],'png');
pause

figure(2);
surf(x,y,u1(:,:,ts),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y');...
    caxis([-1.0,1.0]);...
    text(x(ti),y(ti),10*u1(ti,ti,1),'$u_{1}$','Interpreter','latex','FontSize',fs,'Color','k');...
    set(gca,'xTick',x(1):x(nx)-x(2):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(2):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'u_snapshot'],'png');
pause

qlim1 = min(min(PV_av));
qlim2 = max(max(PV_av));

figure(3);
surf(x,y,PV_av(:,:),'edgecolor','none');view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    caxis([qlim1 qlim2]);
    xlabel('x'); ylabel('y'); title('PV av');...
    set(gca,'xTick',x(1):x(nx)-x(2):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(2):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'PV_av'],'png');
pause

close all

qlim1 = min(min(min(PV1)));
qlim2 = max(max(max(PV1)));

movie = 1;
if movie == 1
    
    for ii=1:nt
        ii
        figure(1)
        surf(PV1(:,:,ii),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet); axis image;...
            xlabel('x'); ylabel('y'); title('PV');...
            caxis([f(ny)/H1,5.0e-7]);
            %caxis([qlim1 qlim2]);...
        
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        %if ii == 1

        %    imwrite(imind,cm,'/home/mike/Documents/GulfStream/GOLD/IMAGES/GOLD_PV.gif', 'Loopcount',inf);

        %else

        %    imwrite(imind,cm,'/home/mike/Documents/GulfStream/GOLD/IMAGES/GOLD_PV.gif','WriteMode','append');
        %end

    end
end

pause 

close all

%%
% 

%loc = '~/cluster/gold3/'; cd(loc);

h = ncread(files(loop_end).name,'h');

h1 = squeeze(h(:,:,1,:));
h2 = squeeze(h(:,:,2,:));
h3 = squeeze(h(:,:,3,:));

h_size = size(h1);
nt = h_size(3);

% h_relax is the relaxation forcing that is applied to h.
h_relax = zeros(h_size);
h_relax_av = zeros(ny,nx);
% h_restore is the constant profile towards which h is being relaxed.
h1_restore = zeros(ny,nx);

for tii = 1:nt
    h1(:,:,tii) = transpose(h1(:,:,tii));
    h2(:,:,tii) = transpose(h2(:,:,tii));
    h3(:,:,tii) = transpose(h3(:,:,tii));
    h0(:,:,tii) = h1(:,:,tii) + h2(:,:,tii) + h3(:,:,tii);
end

% % Relaxation
% %+==========================
% m = 0.1;
% j_shift = 10;
% y_shift = y(j_shift+1) - y(1);
% yf = y(floor(ny/2+1));
% xf = y(floor(ny/2+1));
% L4 = 0.001 * (L/4);
% 
% 
% % relax_const = 1.5e-5;
% % dh = 120.0;
% % for j = 1:ny
% %     for i = 1:nx
% %         if (m * x(i) + yf - 0.5 * L4 < y(j) - y_shift) && (y(j) - y_shift < m * x(i) + yf + 0.5 * L4)
% %             xy = 2.0 * pi * (m * x(i) - y(j) + y_shift + yf) / L4;
% %             sine = sin(xy);
% %             h1_restore(j,i) = 300.0 + dh * sine;
% %             h_relax(j,i,ti) = relax_const * (h1_restore(j,i) - h1(j,i,ti));
% %         %elseif (x(i) < xf) % uncomment for west
% %         elseif (x(i) > xf) % uncomment for east
% %         %else               % uncomment for all
% %            h1_restore(j,i) = 300.0;
% %            h_relax(j,i,ti) = relax_const * (h1_restore(j,i) - h1(j,i,ti));
% %         end
% %     end
% % end
% 
% surf(h1_restore,'edgecolor','none'); view(0,90); colorbar();
% pause

%===========================
h_relax_av = sum(h_relax,3) / nt;

h1lim1 = min(min(min(h1)));
h1lim2 = max(max(max(h1)));
h2lim1 = min(min(min(h2))); 
h2lim2 = max(max(max(h2)));
h3lim1 = min(min(min(h3)));
h3lim2 = max(max(max(h3)));
relaxlim = max(max(max(abs(h_relax_av))));

figure(6)
surf(x,y,h1(:,:,ts),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y');...
    caxis([200.0,400.0]);...
    text(x(ti),y(ti),10*h1(ti,ti,nt),'$h_{1}$','Interpreter','latex','FontSize',fs,'Color','k');...
    set(gca,'xTick',x(1):x(nx)-x(2):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(2):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'h_snapshot'],'png');
pause

figure(7)
surf(x,y,h0(:,:,ts),'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y'); title('ssh');...
    set(gca,'xTick',x(1):x(nx)-x(1):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(1):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'ssh_snapshot'],'png');
pause


for ii=1:nt
    ii
    figure(1)
    surf(h1(:,:,ii),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet); axis image;...
        caxis([h1lim1 h1lim2]);
    mov(ii)=getframe(1);
    im=frame2im(mov(ii));
    [imind,cm] = rgb2ind(im,256);
end

pause 

%%
% Deformation Radius

g12 = 0.01;     % reduced gravity between 1st and 2nd layers
g23 = 0.3 * 0.01;

% new values are 0.01 and 3*0.01

Rd1 = zeros(ny,nx);
Rd2 = zeros(ny,nx);
for i=1:nx
    for j =1:ny
        S1 = f(j)^2 / (g12 * h1(j,i,nt));
        S21 = f(j)^2 / (g12 * h2(j,i,nt));
        S22 = f(j)^2 / (g23 * h2(j,i,nt));
        S3 = f(j)^2 / (g23 * h3(j,i,nt));
        a = S1 + S21 + S22 + S3;
        b = S1 * S22 + S1 * S3 + S21 * S3;
        Rd1(j,i) = 0.001*sqrt((a + sqrt(a^2 - 4 * b)) / (2 * b));
        Rd2(j,i) = 0.001*sqrt((a - sqrt(a^2 - 4 * b)) / (2 * b));
    end
end

figure(7)
surf(x,y,Rd1,'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y'); title('Rd');...
    set(gca,'xTick',x(1):x(nx)-x(1):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(1):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'Rd1_snapshot'],'png');

pause

figure(8)
surf(x,y,Rd2,'edgecolor','none'); view(0,90); shading interp; colorbar(); colormap(jet); axis image;...
    xlabel('x'); ylabel('y'); title('Rd');...
    set(gca,'xTick',x(1):x(nx)-x(1):x(nx)); set(gca,'xTickLabel',{'0','L'});...
    set(gca,'yTick',y(1):y(ny)-y(1):y(ny)); set(gca,'yTickLabel',{'0','L'});...
    saveas(gcf,[outloc,'Rd2_snapshot'],'png');

pause 

close all

%%
% ENERGY

% files = dir('energy__*');
% nf = size(files);
% count = 0;
% %for i = 20:nf
% for i = 1:nf-1
%     disp(i);
%     KEnew = ncread(files(i).name,'KE');
%     KEnew = KEnew(:,:,1,1);
%     nn = size(KEnew,4);    
%     for ti = 1:nn
%         KE1(count+ti+1) =+ sum(sum(KEnew));
%     end
%     count = count + nn;
% end
% nt = size(KE1,3);

%TE = load('energy_d_mass');
%TE = TE(:,2);

%plot(TE); xlabel('days'); ylabel('Total Energy');...
 %   saveas(gcf,['~/Documents/GulfStream/GOLD/Images/','TE'],'png');
    
%pause close all


%%

% Final section can be run independently to find required gprimes to attain
% desired Rossby deformation radii.

g12 = 0.01;     % reduced gravity between 1st and 2nd layers
g23 = 0.3*g12;

%loc = '~/cluster/gold5/';
cd(loc);
files = dir('prog__*');

y = ncread(files(1).name,'yq');
ny = size(y,1);

f0 = 0.44e-4;
beta = 2e-11;
f = f0 + beta * y;

H1 = 300.0;
H2 = 700.0;
H3 = 3000.0;

ff = f(ny-10);
S1 = ff^2 / (g12 * H1);
S21 = ff^2 / (g12 * H2);
S22 = ff^2 / (g23 * H2);
S3 = ff^2 / (g23 * H3);

a = S1 + S21 + S22 + S3;
b = S1 * S22 + S1 * S3 + S21 * S3;

Rd1 = 0.001*sqrt((a + sqrt(a^2 - 4 * b)) / (2 * b))
Rd2 = 0.001*sqrt((a - sqrt(a^2 - 4 * b)) / (2 * b))

%%
close all
y = linspace(0,L,ny);
compU0 = 0;
if compU0 == 1
    U0=zeros(ny,1);
    sigma = 0.02 * L;
    l = L / 2;
    a = max(u_slice) / (exp(l^2 / (2.0 * sigma^2)) - 1.0);
    n2 = ny / 2;
	for j = 1:ny
		U0(j) = a * exp((l^2 - (y(j)-y(n2-3))^2) / (2. * sigma^2)) - a;		% -a ensures U0 is zero on the boundaries
        
    end
    
    plot(y,U0,y,u_slice);
end

%%

% ENERGY

% files = dir('energy__*');
% nf = size(files);
% count = 0;
% %for i = 20:nf
% for i = 1:nf-1
%     disp(i);
%     KEnew = ncread(files(i).name,'KE');
%     KEnew = KEnew(:,:,1,1);
%     nn = size(KEnew,4);    
%     for ti = 1:nn
%         KE1(count+ti+1) =+ sum(sum(KEnew));
%     end
%     count = count + nn;
% end
% nt = size(KE1,3);

TE = load('energy_d_mass');
TE = TE(:,2);

tyears = (1:length(TE)) / 365;

end_ = 40 * 365;
yr10 = 10 * 365 + 1;

E_av = sum(TE(yr10:end_)) / (end_ - yr10)

% Plot full time series
plot(tyears(1:end_),TE(1:end_)/E_av); xlabel('Years'); ylabel('Total Energy / Mass');...
    xlim([0,40]);...
    saveas(gcf,[outloc,'TE'],'png');

pause

% Plot shortened time series
figure('rend','painters','pos',[10 10 700 300])
plot(tyears(yr10:end_),TE(yr10:end_)/E_av); xlabel('Years'); ylabel('Total Energy / Mass');...
    xlim([10,40]);...
    saveas(gcf,[outloc,'TE2'],'png');

    
%pause close all






