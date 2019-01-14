function S = find_length(xx, yy, zz)
% params should have x0, a (slope of x as function of z), y0, and m (slope
% of y line. [x0 y0 a m] x = x0 = a*z , etc
 

%    """ Get the Singular Value Decomposition of the Event. 
%        This returns a line that represents the average behavior of the event, 
%        like PCA """
top_cm = 50; % cm 
rmax =double(24.5); % cm this is the max radius of lux assuming circular shape
drift_to_cm = 1.51/10 ; % cm/us

zz=double(zz);
event = [xx,yy,zz];
event=double(event);
datamean = mean(event);
x0=datamean(1);
y0=datamean(2);
z0=datamean(3);

%  # Do an SVD on the mean-centered data.
v = pca(event - datamean);
   
a=v(1,1);
b=v(2,1);
c=v(3,1);
    
    %params = [x0,y0,z0,a,b,c];
X0= x0 - a./c .* z0;
A=a./c;
Y0=y0 - b./c .* z0;
B=b./c;
 
syms z; % x y
%eqns = [q == x0 + a .* z , y == y0 + m .* z , q.^2 + y.^2 == rmax.^2];
% for some reason matlab acts weird when you have more eqns and a single
% unknown. So I eliminated them 
eqn =X0.^2 + 2.*A.*X0.*z + A.^2*z.^2 + Y0.^2 + 2.*B.*Y0.*z + B.^2.*z.^2 == rmax.^2;
z = solve(eqn, z); % z is the positions where the particle hists the boundary
z =double (z);

for ii=1:length(z)
    if z(ii) > 32000 
        z(ii) = 32000;
    end
    if z(ii) < 0
        z(ii) = 0;
    end
end


if isempty(z)
    z(1) = 0;
    z(2) = 32000;
end

if length(z) == 1
    
    S = 0;
else
    
    x = X0 + A .* z;
    y = Y0 + B .* z;
    drift_us_z = 10 .*  0.001 .* z;
    z_cm = top_cm - drift_us_z .* drift_to_cm;
    
    S = sqrt((x(1) - x(2)).^2 + (y(1) - y(2)).^2 + (z_cm(1) - z_cm(2))^2);
    S=double(S);
end

if S==0
    S=NaN;
end


end