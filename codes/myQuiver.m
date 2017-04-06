function qiu = myQuiver(p2,p1,hsize)

qiu = quiver3(p2(1),p2(2),p2(3),p1(1)-p2(1),p1(2)-p2(2),p1(3)-p2(3),...
    'k','Autoscale','off','Linewidth',1.5,'MaxHeadSize',hsize);
end