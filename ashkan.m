clc
clear
m=100;
n=10;          %n<=m
cor=rand(m,3);
cor_s=cor(randi([1,m],n,1),:);

%% your function
I=1:size(cor,1);
J=zeros(size(cor_s,1),1);
for i=1:size(cor_s,1)   
    d=bsxfun(@minus,cor,cor_s(i,:));
    ind=sqrt(sum(d.^2,2))<1e-6;
    J(i)=I(ind);
end
J  % thisi is what you want
%%
% Test:
cor_s
cor(J,:)

