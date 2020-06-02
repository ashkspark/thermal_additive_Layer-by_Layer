function show(elements3, elements4, coordinates, u)
%figure;
trisurf(elements3, coordinates(:,1), coordinates(:,2), u',...
    'facecolor','interp')
hold on
trisurf(elements4, coordinates(:,1), coordinates(:,2), u',...
    'facecolor', 'interp')
hold off
view(0,90);
%drawnow;
%set(gcf, 'Position',  [100, 100, 500, 400])
title('Solution of the Problem')
end

