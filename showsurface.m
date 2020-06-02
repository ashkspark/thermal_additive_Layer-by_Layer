function showsurface(surface,coordinates,u)
trisurf(surface,coordinates(:,1),coordinates(:,2),...
    coordinates(:,3),u','facecolor','interp')
axis off
colorbar;
%view(0,90)
end

