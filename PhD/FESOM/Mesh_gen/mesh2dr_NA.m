% Channel: regular mesh
% Update: numbering in ijk flavour
  dx=1/4;    % channel u_wide 1/8
  %ddx = 1/8; %in case we want to increase the resolution in a certain part of the mesh
  %dy=dx;
  dy=1/4;   % channel u_wide 1/10
  %ddy = 1/8; 
  lon=0:dx:40;
  lat=30:dy:60; %we may use a higher resolution in the middle of the domain 
  
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
  
  % Set a western coastline and remove all nodes lying on land
 %{
  i=1;
  while i < length(xcoord)
      if ycoord(i) < 60 - 6*xcoord(i) %xcoord(i) < 10 - 0.1*(ycoord(i) - 45).^2
          ycoord(i) = NaN;
          xcoord(i) = NaN;
          i = i + 1;
      else
          i = i + 1;
      end
  end
  %}
  
  tri=[];
  for n=1:nx-1,
      for nn=1:ny-1
          if ~(isnan(xcoord(nodnum(nn,n))) |  isnan(xcoord(nodnum(nn+1,n))) | isnan(xcoord(nodnum(nn,n+1))))
            tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
          end
          if ~(isnan(xcoord(nodnum(nn+1,n))) |  isnan(xcoord(nodnum(nn+1,n+1))) | isnan(xcoord(nodnum(nn,n+1))))
            tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
          end
      end;
  end;    
  
  %trimesh(tri, xcoord', ycoord', zeros(nx*ny,1)); 
  %view(2)
  
  TRI = tri(:);
  
  %remove nodes connected to only one triangle
  
  while min(equaltest(TRI,TRI)) < 2
      
      for i = 1:length(xcoord)
        S = sum(TRI == i);
        if S <= 1 && ~isnan(xcoord(i))
            xcoord(i)
            ycoord(i)
            xcoord(i) = NaN;
            ycoord(i) = NaN;
        end
      end 
      
      tri=[];
      for n=1:nx-1,
          for nn=1:ny-1
              if ~(isnan(xcoord(nodnum(nn,n))) |  isnan(xcoord(nodnum(nn+1,n))) | isnan(xcoord(nodnum(nn,n+1))))
                tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
              end
              if ~(isnan(xcoord(nodnum(nn+1,n))) |  isnan(xcoord(nodnum(nn+1,n+1))) | isnan(xcoord(nodnum(nn,n+1))))
                tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
              end
          end;
      end;    

      TRI = tri(:);
  end
  
  
  
  %tri=delaunay(xcoord, ycoord);
  
  % plot
  figure(1)
  trimesh(tri, xcoord, ycoord, zeros(nx*ny,1));
  xlabel('Longitude') 
  ylabel('Latitude')
  view(2)
  
  % Cyclic reduction:
  % The last ny nodes in xcoord, ycoord are equivalent to 
  % the first ny nodes.
  %ai=find(tri>(nx-1)*ny);
  %tri(ai)=tri(ai)-(nx-1)*ny;
  %xcoord=xcoord(1:(nx-1)*ny);
  %ycoord=ycoord(1:(nx-1)*ny);
 
  % Cyclic reduction means that the last column has to be removed 
  % in xnum, ynum 
  %xnum=xnum(:, 1:nx-1);
  %ynum=ynum(:, 1:nx-1);
  
  %xnum=reshape(xnum,[ny*(nx-1),1]);
  %ynum=reshape(ynum,[ny*(nx-1),1]);
  
 
  n2d=(nx)*ny;
  nodes=zeros([4, n2d]);
  
  
  nodes(1,:)=1:n2d;
  nodes(4,:)=zeros(size(ycoord));
  % Set indices to 1 on vertical walls
  %ai=find(ycoord==min(lat));
  %nodes(4,ai)=1;
  %ai=find(ycoord==max(lat));
  %nodes(4,ai)=1;
  %ai=find(xcoord==min(lon));
  %nodes(4,ai)=1;
  %ai=find(xcoord==max(lon));
  %nodes(4,ai)=1;
  
  
  %Set indices to 1 on boundary:
  
  TRI = tri(:);
  for i=1:size(xcoord,2)
      %How many times does node i appear in TRI? That is, how many
      %triangles is node i attached to?
      ai = sum(TRI == i);
      if ai < 6
         nodes(4,i) = 1 ;
      end
  end
  
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
  
  figure(2)
  plot3(nodes(2,:),nodes(3,:),nodes(4,:),'.');
  
  %remove the two edge nodes 1 and nx*ny with respective elements with only one neighbor
  %tri=tri(2:end-1,:);
  %nodes=nodes(:, 2:end-1);
  
  %Relabel nodes so that they take numbers 1,2,3,...
  newnodnum = int64(isnan(xcoord));
  newlen = length(xcoord) - sum(newnodnum);
  newnodnum = 1-newnodnum;
  ai = find(newnodnum == 1);
  i=1;
  for a=ai
      newnodnum(a) = i;
      i = i+1;
  end
  
  n2d = newlen;
  
  %Remove NaN nodes from xcoords and ycoords
  
  i=1;
  while i <= length(xcoord)
      if isnan(xcoord(i))
          xcoord(i) = [];
          ycoord(i) = [];
          nodes(:,i) = [];
      else
          i = i+1;
      end
  end
  
  %Convert triangles to new node numbering
  
  tri = newnodnum(tri);
  size(nodes)
  nodes(1,:) = newnodnum(nodes(1,:));
 
  
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
      %{
      shelflength = 2.5;
      shelfdepth = -100;
      slopeend = 10;
      gradient = (-shelfdepth - 1600)/(slopeend - shelflength);
      c = shelfdepth - gradient*shelflength;
      depth= shelfdepth*(xcoord < shelflength) - 1600*(xcoord>slopeend) + (xcoord>=shelflength).*(xcoord<=slopeend).*(gradient*xcoord + c) ;
      %}
      
      fid=fopen('aux3d.out', 'w');
      fprintf(fid,'%g\n', nl);
      fprintf(fid,'%g\n', zbar);
      fprintf(fid,'%7.1f\n', depth);
      fclose(fid);
     
 %fid=fopen('xynum.out','w');
 %       fprintf(fid,'%8i %8i\n',[xnum, ynum]);
 %       fclose(fid);
  
 

