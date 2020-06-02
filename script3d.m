% points(x(:),y(:),z(:));
% Unrecognized function or variable 'points'.
clear all; clc;
x = linspace(0,4,5);
y = linspace(0,4,5);
z = linspace(0,5,6);
LN = 5;
[x, y, z] = meshgrid(x,y,z);
points = [x(:),y(:),z(:)];
num = size(points,1);
DT = delaunayTriangulation(points);
Tnum = size(DT,1);
trimesh(DT.ConnectivityList,x,y,z)
u = 1:num;
Tu = 1:Tnum;
pointsc = [u',points,LN*ones(num,1)];  
DTc = [Tu',DT.ConnectivityList,LN*ones(Tnum,1)]; 
idx1 = points(:,3) == 0; %points with z=0 in order to find the dirichlet surface
us=u';
dirich = us(idx1,:);
[F,P] = freeBoundary(DT);
J = zeros(size(P,1),1);
%%
for i=1:size(P,1)   
    d = bsxfun(@minus,DT.Points,P(i,:));
    ind=sqrt(sum(d.^2,2))<1e-6;
    J(i)=u(ind);
   % F(F==i) = u(ind);
end

for kk =1:size(F(:,1),1)
    F(kk,1) = J(F(kk,1));  
    F(kk,2) = J(F(kk,2)); 
    F(kk,3) = J(F(kk,3)); 
end

TF = size(F,1);
inXd = zeros(TF,1);
%-------------------------------
for ss=1:TF
    inXd(ss) = ismember(F(ss,1),dirich) && ismember(F(ss,2),dirich) && ismember(F(ss,3),dirich);
end
%-------------------------------
inXd = boolean(inXd);
Fneum = F(~inXd,:);
Fdirich = F(inXd,:);
nne = size(Fneum,1);
nnd = size(Fdirich,1);
une = 1:nne;
und = 1:nnd;
Fneum = [une',Fneum,LN*ones(nne,1)];
Fdirich = [und',Fdirich,LN*ones(nnd,1)];
%uF  = unique(F);