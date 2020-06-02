function hottemp = Th(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
hottemp = 50*ones(size(x,1),1);
idd = (x(:,1)-2).^2+(x(:,2)-2).^2 <= 2.0; 
hottemp(idd,:) = 100;

%ones(size(x,1),1)
end

