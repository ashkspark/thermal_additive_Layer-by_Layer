%FEM3D_HEAT finite element method for three-dimensional additive manufacturing heat equation.
%Copyright Ashkan Golgoon, PhD 
clear all; clc;
%%
%Loading the Data
load coordinates_heat.dat;
load elements3_heat.dat;
load neumann_heat.dat;
load dirichlet_heat.dat;
%%
T = 1; dL = 0.20; NL = T/dL; dt = 0.001; %time needed to create each layer (could be variable among layers)
N = dL/dt;  %number of layers
%---------------------------
%Th = 100*ones(1,NL); %hot layer temparatures
%Th(2) = 850;
fra = 4;
%---------------------------
num = zeros(1,NL);
nt  = zeros(1,NL);
ne  = zeros(1,NL);
nd  = zeros(1,NL);
%%
coordinates_heat(:,1)=[]; %empty out the first column
elements3_heat(:,1)=[];
neumann_heat(:,1) = [];
dirichlet_heat(:,1) = [];
for i = 1:NL
%---------------------------------    
%size(coordinates_heat(:,1))
indxc = coordinates_heat(:,4) == i;
cor{i} = coordinates_heat(indxc,:);
cor{i}(:,4) = []; %dropping the layer index
num(i)= size(coordinates_heat(indxc,:),1);
%---------------------
%Initialization
A{i} = sparse(num(i),num(i));
B{i} = sparse(num(i),num(i));
b{i} = sparse(num(i),1);
%num = size(coordinates_heat,1);
%---------------------------------
indxe = elements3_heat(:,5) == i; 
elem{i} = elements3_heat(indxe,:);
elem{i}(:,5) = [];
nt(i) = size(elements3_heat(indxe,:),1);
%---------------------------------
indxn = neumann_heat(:,4) == i; 
neum{i} = neumann_heat(indxn,:);
neum{i}(:,4) = [];
ne(i) = size(neumann_heat(indxn,:),1);
%---------------------------------
indxd = dirichlet_heat(:,4) == i;  
dirich{i} = dirichlet_heat(indxd,:);
dirich{i}(:,4) = [];
nd(i) = size(dirichlet_heat(indxd,:),1);
%---------------------------------
FreeNodes{i} = setdiff(1:num(i),unique(dirich{i})); %Finding the free nodes
%Initializing the rhs and the stiffness matrix
%A = sparse(num,num,NL); %all zero sparse matrix, stiffness matrix
%B = sparse(num,num,NL); %all zero sparese matrix, mass matrix
%Assembly
%assembling the triangular mesh
for  j = 1:nt(i)
 %stiffness matrix   
    A{i}(elem{i}(j,:),elem{i}(j,:)) = A{i}(elem{i}(j,:), ...
     elem{i}(j,:)) + stima3(cor{i}(elem{i}(j,:),:));
 %mass matrix
    B{i}(elem{i}(j,:),elem{i}(j,:)) = B{i}(elem{i}(j,:), ...
      elem{i}(j,:)) + abs(det([1,1,1,1;cor{i}(elem{i}(j,:),:)']))...
        *[2,1,1,1;1,2,1,1;1,1,2,1;1,1,1,2]/120;
end
U{i} = zeros(num(i),N+1);
end

%%
for i = 1:NL
%---------------------------------------    
% time steps
%---------------------------------------
for n = 2:N+1
b{i} = sparse(num(i),1);
% Volume Forces
for j = 1:nt(i)
b{i}(elem{i}(j,:)) = b{i}(elem{i}(j,:)) + ...
abs(det([1,1,1,1; cor{i}(elem{i}(j,:),:)'])) * ...
dt*fheat(sum(cor{i}(elem{i}(j,:),:))/4,(i-1)*dL+(n-1)*dt)/24;
end
% Neumann conditions
for j = 1:ne(i)
b{i}(neum{i}(j,:)) = b{i}(neum{i}(j,:)) + ...
norm(cross(cor{i}(neum{i}(j,3),:)-...
cor{i}(neum{i}(j,1),:),cor{i}(neum{i}(j,2),:)-cor{i}(neum{i}(j,1),:)))*...
dt*gheat(sum(cor{i}(neum{i}(j,:),:))/3,(i-1)*dL+(n-1)*dt)/6;
end
% previous timestep
b{i} = b{i} + B{i} * U{i}(:,n-1);
% Dirichlet conditions
u = sparse(num(i),1);
u(unique(dirich{i})) = u_d_heat(cor{i}(unique(dirich{i}),:),(i-1)*dL+(n-1)*dt);
b{i} = b{i} - (dt * A{i} + B{i}) * u;
% Computation of the solution
u(FreeNodes{i}) = (dt*A{i}(FreeNodes{i},FreeNodes{i})+ ...
B{i}(FreeNodes{i},FreeNodes{i}))\b{i}(FreeNodes{i});
U{i}(:,n) = u;
end
%---------------------------------------
% Initial Condition: passing nodal temparature from previous to the next
% step by setting the initial condition
if i <= NL-1
U{i+1}(1:num(i),1) = U{i}(:,n);
U{i+1}(num(i)+1:num(i+1),1) = Th(cor{i+1}(num(i)+1:num(i+1),:));
% U{i+1}(num(i)+3:num(i+1),1) = 0;
end
%---------------------------------------
end
 %%
% graphic representation
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
%vidfile.FrameRate = fra;
open(vidfile);
for i = 1:NL
for kk = 1:fra:N+1
showsurface([dirich{i};neum{i}],cor{i},full(U{i}(:,kk)));
movieVector{i}(kk) = getframe(gcf);
end
end
for i = 1:NL
writeVideo(vidfile,movieVector{i}(1:fra:N+1));
end
close(vidfile);
%for i = 1:NL
% for kk = 1:N+1
% show(elem{4},[],cor{4},full(U{i}(:,kk)));
% movieVector(kk) = getframe;
% end
%end


