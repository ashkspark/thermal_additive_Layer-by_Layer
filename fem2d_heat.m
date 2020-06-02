%FEM2D_HEAT finite element method for two-dimensional heat equation.
clear all; clc;
%Initialisation
load coordinates_heat.dat; coordinates_heat(:,1)=[]; %empty out the first column
num = size(coordinates_heat,1);
load elements3_heat.dat; elements3_heat(:,1)=[];
nt = size(elements3_heat,1);
load neumann_heat.dat; neumann_heat(:,1) = []; 
ne = size(neumann_heat,1);
load dirichlet_heat.dat; dirichlet_heat(:,1) = [];
nd = size(dirichlet_heat,1);
%%
FreeNodes = setdiff(1:num,unique(dirichlet_heat)); %Finding the free nodes
%Initializing the rhs and the stiffness matrix
A = sparse(num,num); %all zero sparse matrix, stiffness matrix
B = sparse(num,num); %all zero sparese matrix, mass matrix
T = 1; dt = 0.01; N = T/dt;
U = zeros(num,N+1);
%%
%Assembly
%assembling the triangular mesh
for  j = 1:nt
 %stiffness matrix   
    A(elements3_heat(j,:),elements3_heat(j,:)) = A(elements3_heat(j,:), ...
     elements3_heat(j,:)) + stima3(coordinates_heat(elements3_heat(j,:),:));
 %mass matrix
    B(elements3_heat(j,:),elements3_heat(j,:)) = B(elements3_heat(j,:), ...
      elements3_heat(j,:)) + det([1,1,1;coordinates_heat(elements3_heat(j,:),:)'])...
        *[2,1,1;1,2,1;1,1,2]/24;
end
%%
% Initial Condition
U(:,1) = zeros(num,1);
%%
% time steps
for n = 2:N+1
b = sparse(num,1);
% Volume Forces
for j = 1:nt
b(elements3_heat(j,:)) = b(elements3_heat(j,:)) + ...
det([1,1,1; coordinates_heat(elements3_heat(j,:),:)']) * ...
dt*fheat(sum(coordinates_heat(elements3_heat(j,:),:))/3,n*dt)/6;
end
% Neumann conditions
for j = 1:ne
b(neumann_heat(j,:)) = b(neumann_heat(j,:)) + ...
norm(coordinates_heat(neumann_heat(j,1),:)-coordinates_heat(neumann_heat(j,2),:))*...
dt*gheat(sum(coordinates_heat(neumann_heat(j,:),:))/2,n*dt)/2;
end
% previous timestep
b = b + B * U(:,n-1);
% Dirichlet conditions
u = sparse(num,1);
u(unique(dirichlet_heat)) = u_d_heat(coordinates_heat(unique(dirichlet_heat),:),n*dt);
b = b - (dt * A + B) * u;
% Computation of the solution
u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
B(FreeNodes,FreeNodes))\b(FreeNodes);
U(:,n) = u;
end
 %%
% graphic representation
for kk = 1:N+1
show(elements3_heat,[],coordinates_heat,full(U(:,kk)));
movieVector(kk) = getframe;
end





