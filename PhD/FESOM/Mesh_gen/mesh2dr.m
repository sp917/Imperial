% Channel: regular mesh
% Update: numbering in ijk flavour
  dx=1/12;    % channel u_wide 1/8
  %dy=dx;
  dy=1/14;   % channel u_wide 1/10
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
   
  xcoord=reshape(xcoord,[1,nx*ny]);
  ycoord=reshape(ycoord,[1,nx*ny]); 
  
  tri=[];
  for n=1:nx-1,
      for nn=1:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
      tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
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
  nodes(4,:)=zeros(size(ycoord));
  % Set indices to 1 on vertical walls
  ai=find(ycoord==min(lat));
  nodes(4,ai)=1;
  ai=find(ycoord==max(lat));
  nodes(4,ai)=1;

  % Mesh distortion:
  amp=0.4;
  ampx=amp*dx;
  ampy=amp*dy;
  if 1<0,
     ai=find(nodes(4,:)==0);
     xcoord(ai)=xcoord(ai)+(ampx*(-0.5+rand([1, length(ai)])));
     ycoord(ai)=ycoord(ai)+(ampy*(-0.5+rand([1, length(ai)]))); 
  end;
 
  nodes(2,:)=xcoord;
  nodes(3,:)=ycoord;
  
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
 
 % Define levels:
      zbar=[0 10 22 35 49 63 79 100 150 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600]; 
      %zbar=[0:100:1500];
      nl=length(zbar);
      % define depths:
      %depth=-500-3500*sin((ycoord-30)*pi/15.0);
      depth=-1600*ones(size(xcoord));
      %depth=-1500*ones(size(xcoord));
      
      fid=fopen('depth.out', 'w');
      fprintf(fid,'%g\n', nl);
      fprintf(fid,'%g\n', zbar);
      fprintf(fid,'%7.1f\n', depth);
      fclose(fid);
     
 %fid=fopen('xynum.out','w');
 %       fprintf(fid,'%8i %8i\n',[xnum, ynum]);
 %       fclose(fid);
  
