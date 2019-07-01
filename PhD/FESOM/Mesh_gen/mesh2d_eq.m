% Channel: regular mesh, equilateral triangles
% Update: numbering in ijk flavour
  dx=1/6;
  dy=dx*cos(37.5*pi/180);    
  dy=dy*sqrt(3)/2;
  lon=0:dx:20;
  lat=30:dy:45;
  
  nx=length(lon);
  ny=length(lat);
  nodnum=1:nx*ny;
  nodnum=reshape(nodnum,[ny, nx]);
  xcoord=zeros([ny, nx]);
  xnum=zeros([ny, nx]);
  ycoord=xcoord;
  ynum=zeros([ny, nx]);
  for n=1:nx,
  ycoord(:,n)=lat';
  ynum(:,n)=(1:ny)';
  end;
   
  for n=1:ny,
  xcoord(n,:)=lon;
  xnum(n,:)=(1:nx);
  end;
  for n=2:2:ny,
  xcoord(n,:)=xcoord(n,:)+0.5*dx;
  end;
  xcoord=reshape(xcoord,[1,nx*ny]);
  ycoord=reshape(ycoord,[1,nx*ny]); 
  
  tri=[];
  for n=1:nx-1,
      for nn=1:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
      tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
      end;
      for nn=2:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn+1,n+1)]];
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
      end;
  
  
  end;    
  %tri=delaunay(xcoord, ycoord);
  
  % plot
  trimesh(tri, xcoord', ycoord', zeros(nx*ny,1)); 
  view(2)
  
 
  % Cyclic reduction:
  % The last ny nodes in xcoord, ycoord are equivalent to 
  % the first ny nodes.
  ai=find(tri>(nx-1)*ny);
  tri(ai)=tri(ai)-(nx-1)*ny;
  xcoord=xcoord(1:(nx-1)*ny);
  ycoord=ycoord(1:(nx-1)*ny);
 
  % Cyclic reduction means that the last column has to be removed 
  % in xnum, ynum 
  xnum=xnum(:, 1:nx-1);
  ynum=ynum(:, 1:nx-1);
  xnum=reshape(xnum,[ny*(nx-1),1]);
  ynum=reshape(ynum,[ny*(nx-1),1]);
  
  n2d=(nx-1)*ny;
  nodes=zeros([4, n2d]);
  nodes(1,:)=1:n2d;
  nodes(2,:)=xcoord;
  nodes(3,:)=ycoord;
  nodes(4,:)=zeros(size(ycoord));
  % Set indices to 1 on vertical walls
  ai=find(ycoord==min(lat));
  nodes(4,ai)=1;
  ai=find(ycoord==max(lat));
  nodes(4,ai)=1;

  % Define levels:
      levels=[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600];  
      %levels=[0:100:1500];
      nl=length(levels);
      % define depths:
      %depth=-500-3500*sin((ycoord-30)*pi/15.0);
      dd=-1600*ones(size(xcoord));
      %dd=-1500*ones(size(xcoord));
      %dd=-1600+500*exp(-((xcoord-10).^2+(ycoord-37.5).^2)/9);
      elem=tri';
      e2d=length(tri(:,1));
      nlevels=zeros([1,e2d]);
      lmid=[0.5*(levels(1:end-1)+levels(2:end)),levels(end)];
     for nt=1:e2d,
      ai=find(lmid+sum(dd(elem(:,nt)))/3.0>=0);   
      nlevels(nt)=ai(1);
     end;
ai=find(nlevels<3); nlevels(ai)=3;
%figure(4)
%patch(xcoord(elem), ycoord(elem), nlevels);
      
      
%mesh_circle_aux
%nlcopy(1)=24;
%nlcopy(24618)=24;  

%dhb=sum(dd(elem))/3.0+levels(nlcopy);
%ai=find(dhb>=0);  % model depth is too large
%dhb(ai)=min(dhb(ai),(levels(nlcopy(ai))-levels(nlcopy(ai)-1))/2);

%ai=find(dhb<0);   % model depth is too small
%dhb(ai)=max(dhb(ai),-(levels(nlcopy(ai))-levels(nlcopy(ai)-1))/2);
%dhb=-dhb;   % This has to be added to dz to have actual thickness

  % Output 2D mesh 
  fid = fopen('nod2d.out','w');
        fprintf(fid,'%8i \n',n2d);
        fprintf(fid,'%8i %8.4f %8.4f %8i\n',nodes); 
        fclose(fid);
 clear nodes       
        
 fid=fopen('elem2d.out','w');
        fprintf(fid,'%8i \n', length(tri(:,1)));
        fprintf(fid,'%8i %8i %8i\n',tri');
        fclose(fid);
 
 
      fid=fopen('depth.out', 'w');
      fprintf(fid,'%g\n', nl);
      fprintf(fid,'%g\n', levels);
      fprintf(fid,'%7.1f\n', dd);
      %fprintf(fid,'%8i\n',nlcopy); 
      %fprintf(fid,'%7.1f\n', dhb);
       
fclose(fid);


