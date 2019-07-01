dx = 1/4;
dy = 1/4;

lon = 0:dx:40;
lat = 30:dy:60;

nx = length(lon);
ny = length(lat);

xcoord = zeros(ny,nx);
ycoord = zeros(ny,nx);

for n=1:ny
    for nn=1:nx
        xcoord(n,nn) = lon(nn);
        ycoord(n,nn) = lat(n);
    end
end

%Define a western coastline and remove all nodes on land

xcoord=reshape(xcoord,[1,nx*ny]);
ycoord=reshape(ycoord,[1,nx*ny]); 

% Set a western coastline and remove all nodes lying on land
for i = 1:length(xcoord)
   if xcoord(i) < 10 - 0.1*(ycoord(i) - 45).^2
      ycoord(i) = NaN;
      xcoord(i) = NaN;
   end
end

size(xcoord)
size(ycoord)

xcoord = rmmissing(xcoord);
ycoord = rmmissing(ycoord);

size(xcoord)
size(ycoord)








