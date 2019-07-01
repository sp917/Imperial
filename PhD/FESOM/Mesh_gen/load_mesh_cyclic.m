
% Load the mesh:
meshdir='./';
cyclic=0;
cyclic_length=20;
%mesh_loaded=1;
%if exist('mesh_loaded')==0,
 mesh_loaded=1;  

 fid=fopen([meshdir,'nod2d.out']);
 n2d=fscanf(fid, '%g',1);
 nodes=fscanf(fid, '%g',[4, n2d]);
 xcoord=nodes(2,:);
 ycoord=nodes(3,:);
 nodind=nodes(4,:);
 clear nodes
 fclose(fid);
 fid=fopen([meshdir,'elem2d.out']);
 e2d=fscanf(fid, '%g',1);
 elem=fscanf(fid, '%g',[3, e2d]);
 fclose(fid);
 
 %fid=fopen([meshdir,'depth.out']);
 fid=fopen([meshdir,'aux3d.out']);
 nl=fscanf(fid, '%g',1);
 zbar=fscanf(fid, '%g',[1,nl]);
 depth=fscanf(fid, '%g',[1,n2d]);
 nlevels=fscanf(fid, '%g',[1,e2d]);
 dhb=fscanf(fid, '%g',[1,e2d]);
 fclose(fid);
 %ZZ=(zbar(1:nl-1)+zbar(2:nl))/2;  % midlevels
 xc=xcoord(elem);
 yc=ycoord(elem);
    
 if cyclic,
    if 1>2,
    xmin=min(xc);
    x1=xc(1,:)-xmin;
    x2=xc(2,:)-xmin; 
    x3=xc(3,:)-xmin;
    
    ai=find(x1>cyclic_length/2);
    xc(1,ai)=xc(1,ai)-cyclic_length;
    ai=find(x2>cyclic_length/2);
    xc(2,ai)=xc(2,ai)-cyclic_length;
    ai=find(x3>cyclic_length/2);
    xc(3,ai)=xc(3,ai)-cyclic_length;
    end;
    xmax=max(xc);
    x1=xc(1,:)-xmax;
    x2=xc(2,:)-xmax; 
    x3=xc(3,:)-xmax;
    
    ai=find(x1<-cyclic_length/2);
    xc(1,ai)=xc(1,ai)+cyclic_length;
    ai=find(x2<-cyclic_length/2);
    xc(2,ai)=xc(2,ai)+cyclic_length;
    ai=find(x3<-cyclic_length/2);
    xc(3,ai)=xc(3,ai)+cyclic_length;
 end;  
% end;
%figure(2); patch(xc,yc,'g') 
figure; patch(xc,yc,'g') 