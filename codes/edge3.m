addpath('D:\TDV\tdv_toolbox')
addpath C:\Users\kubin\Desktop\GitFolder\doc_thesis\codes
fig_path = 'C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\';

% definice bbx
bbx = [5 -4 15 4];
bbx2 = [-7 -4 15 25];
% bbx2 = [2 -4 24 25];

p1  = [bbx(1);bbx(4);1];
p2  = [bbx(3);bbx(4);1];
p3  = [bbx(1);bbx(2);1];
p4  = [bbx(3);bbx(2);1];
% bbx_points = [p1 p2 p4 p3 p1]; 

% bbx lines 
top_line= cross(p1,p2);
bottom_line = cross(p3,p4);
l_line  = cross(p1,p3);
r_line  = cross(p2,p4);

p1  = [bbx2(1);bbx2(4);1];
p2  = [bbx2(3);bbx2(4);1];
p3  = [bbx2(1);bbx2(2);1];
p4  = [bbx2(3);bbx2(2);1];
% bbx_points = [p1 p2 p4 p3 p1]; 

% bbx lines 
top_line2    = cross(p1,p2);
bottom_line2 = cross(p3,p4);
l_line2      = cross(p1,p3);
r_line2      = cross(p2,p4);



figure
% plot(bbx_points(1,:),bbx_points(2,:),'b')
 hold on
 axis equal



%----------------------------%
%----------------------------%
%----------------------------%
addpath('D:\TDV\tdv_toolbox')

fi = 7/10*pi; % pi/2 - cca pi
radius_l = 3;
samples = 30;
vector = [0 1 0]'; % rovna faseta


% volba indexu lomu zmeni pripad lomu na hrane 
no = 1.0;
ni = 1.5;
fromGlass = ni > no;

if fromGlass
    vi = [0 1 0]';
else
    vi = [0 -1 0]';
end



% prunik faset
xmax = 11;
ymin = 1;


ref_fi = pi - fi;
intersPoint = [xmax ymin]';

l1 = [bbx2(1) xmax-radius_l; % faseta 1
      ymin  ymin];
l1_line      = cross([bbx2(1);ymin;1],[xmax-radius_l;ymin;1]);
dx = radius_l * cos(ref_fi);
dy = -radius_l * sin(ref_fi);
facet2_line = cross(e2p(intersPoint), e2p([xmax+dx;ymin+dy]));

l1_points = [cross(facet2_line,bottom_line) cross(facet2_line,r_line)];
f2_max = p2e(l1_points);
out =  logical(isinbbx(f2_max, bbx));
f2_max = f2_max(:,find(out,1,'first'));


l2 = [xmax+dx f2_max(1);
      ymin+dy f2_max(2)];
  
hold on
plot(l1(1,:),l1(2,:),'k','Linewidth',1.5);h = plot(l2(1,:),l2(2,:),'k','Linewidth',1.5);
legsource = h;
%----------------------------%
%----------------------------%
%----------------------------%

r = radius_l/tan(ref_fi/2);

angles = 0:ref_fi/(samples-1):ref_fi;

y = r*cos(angles);
x = r*sin(angles);
radiusX = xmax-radius_l +  x;
radiusY = ymin-r + y;

h = plot(radiusX,radiusY,'m','Linewidth',1.5);
legsource=[legsource h];
plot([bbx2(1);bbx2(3)],[bbx2(4);bbx2(4)],'Color',[0.6 0.6 0.6],'Linewidth',2.5);
axis off
set(gcf,'Color','w')
axis equal

% body hran
xEdge =xmax-radius_l +  x(1:end-1) + diff(x)/2;
yEdge =ymin-r +  y(1:end-1) + diff(y)/2;

fi = atan2(diff(y),diff(x));

n_f = zeros(3,length(fi));
for r = 1:length(fi)
    if fromGlass
        n_f(:,r) = Rz(fi(r) +pi)*vector;
    else
        n_f(:,r) = Rz(fi(r))*vector;
    end
end

vo = ray_reflection(n_f(:,1),vi);

x = xmax-radius_l-0.5 : 0.15: xmax-radius_l;
y = ymin; 
for r = 1:length(x)
     if fromGlass
        quiver(x(r),bbx2(2),0,ymin-bbx2(2),'b','AutoScale','off');
     else
        quiver(x(r),bbx2(4),0,ymin-bbx2(4),'b','AutoScale','off');
        quiver(x(r),ymin,0,bbx2(4)-ymin,'r','AutoScale','off');
     end
     if fromGlass
        quiver(x(r),ymin,0,bbx2(4)-ymin,'g','AutoScale','off');
     else
        quiver(x(r),ymin,0,bbx2(2)-ymin,'g','AutoScale','off');
     end
end


edge_points  = [xEdge;yEdge;ones(1,length(yEdge))];

vect_points = edge_points - 20*repmat(vi,1,length(yEdge));
rayLine   = zeros(3,length(yEdge));
rayPoint = zeros(3,length(yEdge));
mask = 1:length(yEdge);
maskOk = true(1,length(yEdge));
for w = mask
    rayLine(:,w) = cross(vect_points(:,w),edge_points(:,w));
    rayLine(:,w) = rayLine(:,w)/rayLine(3,w);
    if (vi'*n_f(:,w)/(norm(vi)*norm(n_f(:,w)))) < 0
    
        if fromGlass
            l1_points = [cross(rayLine(:,w),bottom_line) cross(rayLine(:,w),l_line)];
        else
            l1_points = [cross(rayLine(:,w),top_line2) cross(rayLine(:,w),r_line)];
        end
        l1_pointsE = p2e(l1_points);
        out =  logical(isinbbx(l1_pointsE, bbx2));
        rayPoint(:,w) = l1_points(:,find(out,1,'first'));
        rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
    %     plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-b')
        myq =  quiver(rayPoint(1,w),rayPoint(2,w),xEdge(w)-rayPoint(1,w),yEdge(w)-rayPoint(2,w),'b','AutoScale','off');
        if length(legsource) < 3
            legsource = [legsource myq];
        end
    else
        maskOk(w) = false;
    end
end
mask = mask(maskOk);

rayLine   = zeros(3,length(mask));
rayPoint = zeros(3,length(mask));
fio   = zeros(1,length(mask)); % odraz
fio1  = zeros(1,length(mask)); % lom
fio2  = zeros(1,length(mask)); 
fio12 = zeros(1,length(mask)); 

cosAlfa = zeros(1,length(mask)); % uhel odrazu
cosBeta = zeros(1,length(mask)); % uhel lomu
for w = mask
    vo = real(ray_reflection(n_f(:,w),vi));
    vectOut = edge_points(:,w) - 20*vo;
    fio(w)  = abs(atan2(vo(1),vo(2))*180/pi);
    fio1(w) = abs(atan2(-vo(1),-vo(2))*180/pi);
    cosAlfa(w)   =  vo'*n_f(:,w)/(norm(vo)*norm(n_f(:,w)));
    
    rayLine(:,w) = cross(vectOut,edge_points(:,w));
    rayLine(:,w) = rayLine(:,w)/rayLine(3,w);
    
    if fromGlass
        l1_points = [cross(rayLine(:,w),bottom_line) cross(rayLine(:,w),l_line2) cross(rayLine(:,w),l1_line)];
        bbx3 = [bbx2(1) bbx2(2) xmax-radius_l ymin];
    else
        l1_points = [cross(rayLine(:,w),top_line2) cross(rayLine(:,w),r_line2)  cross(rayLine(:,w),bottom_line)];
        bbx3 = bbx2;
    end
    l1_pointsE = p2e(l1_points);
    out =  logical(isinbbx(l1_pointsE, bbx3));
    
    rayPoint(:,w) = l1_points(:,find(out,1,'first'));
    rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
    
    if find(out,1,'first') == 3 && fromGlass
         vo = real(ray_reflection([0;0;-1],rayLine(:,w)));
        vectOut2 = rayPoint(:,w) - 20*vo;
 
        rayLine2 = cross(vectOut2,rayPoint(:,w));
        rayLine2 = rayLine2/rayLine2(3);
    
        l1_points = [cross(rayLine2,bottom_line) cross(rayLine2,l_line2) cross(rayLine2,l1_line)];

        l1_pointsE = p2e(l1_points);
        bbx3 = [bbx2(1) bbx2(2) xmax-radius_l ymin];
        out =  logical(isinbbx(l1_pointsE, bbx3));

        rayPoint2 = l1_points(:,find(out,1,'first'));
        rayPoint2 = rayPoint2/rayPoint2(3);
        myq =  quiver(rayPoint(1,w),rayPoint(2,w),-rayPoint(1,w)+rayPoint2(1),-rayPoint(2,w)+rayPoint2(2),'r','AutoScale','off');
    end
    
    plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-r')
    myq =  quiver(xEdge(w),yEdge(w),-xEdge(w)+rayPoint(1,w),-yEdge(w)+rayPoint(2,w),'r','AutoScale','off');
    
    if length(legsource) < 4
        legsource = [legsource myq];
    end
    
    vo = real(ray_refraction(n_f(:,w),vi,ni,no));
    cosBeta(w)   =  vo'*-n_f(:,w)/(norm(vo)*norm(-n_f(:,w)));
    
        
        fio2(w)  = abs(atan2(vo(1),vo(2))*180/pi);
        fio12(w) = abs(atan2(-vo(1),-vo(2))*180/pi);

        vectOut = edge_points(:,w) + 20*vo;

        rayLine(:,w) = cross(vectOut,edge_points(:,w));
        rayLine(:,w) = rayLine(:,w)/rayLine(3,w);

        if fromGlass
           l1_points = [cross(rayLine(:,w),top_line2) cross(rayLine(:,w),r_line) cross(rayLine(:,w),l_line2)];
        else
           l1_points = [cross(rayLine(:,w),bottom_line) cross(rayLine(:,w),r_line)];
        end
        l1_pointsE = p2e(l1_points);
        out =  logical(isinbbx(l1_pointsE, bbx2));
        rayPoint(:,w) = l1_points(:,find(out,1,'first'));
        rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
    %     plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-r')
    if abs(cosBeta(w))-1e-10 > 0 % totalni obraz
        myq = quiver(xEdge(w),yEdge(w),-xEdge(w)+rayPoint(1,w),-yEdge(w)+rayPoint(2,w),'g','AutoScale','off');
        if length(legsource) < 5
            legsource = [legsource myq];
        end
    end
end

if ~fromGlass
    h = legend(legsource,'fasety','hrana',sprintf('vstupující\npaprsek'),'odraz','lom');
    h.FontSize = 12;
    h.Position = [0.1835 0.4763 0.1597 0.2218];
end
text(5.6,1.6,'1','Interpreter','Latex','FontSize',16)
text(14.1,-2.9,'2','Interpreter','Latex','FontSize',16)


set(gcf,'Position',[403 108 774 558])
if fromGlass
%     text(13,-4,'$n_1 = 1.5$','Interpreter','Latex')
    text(-4,-2,'Sklo','FontSize',14) 
    text(-4, 3,'Vzduch','FontSize',14) 
    text(9.2, 26,'Stínítko','FontSize',14) 
    text(11.2, 22.5,'Stopa','FontSize',12) 
    text(-6.2, 22.5,'Ocásek','FontSize',12)
    quiver(-3.5, 23, 3,1.5,'Autoscale','off','Color','k','linewidth',1.5)
    quiver(10.9, 22.5, -2.2,2.0,'Autoscale','off','Color','k','linewidth',1.5)
    
    plotEllipse(xmax-radius_l-0.2,bbx2(4),0.5,0.5,0,[0.3 0.3 0.3],1,0,1.5)
    p = plotEllipse(0.5,bbx2(4),7,0.5,0,[0.3 0.3 0.3],1,0,1);
    set(p,'LineStyle','--')
    
    export_fig([fig_path 'edgeOutFar.pdf'])
else
    text(2,-2,'Sklo','FontSize',14) 
    text(2, 5,'Vzduch','FontSize',14) 
    text(9.2, 26,'Stínítko','FontSize',14) 
    text(2.2, 22.5,'Stopa','FontSize',12) 
    text(16.5, 22.5,'Ocásek','FontSize',12)
    quiver(17.5, 23, -3,1.5,'Autoscale','off','Color','k','linewidth',1.5)
    quiver(5.2, 22.5, 2.2,2.0,'Autoscale','off','Color','k','linewidth',1.5)
    
    plotEllipse(xmax-radius_l-0.2,bbx2(4),0.5,0.5,0,[0.3 0.3 0.3],1,0,1.5)
    p = plotEllipse(0.5+15.5,bbx2(4),8,0.5,0,[0.3 0.3 0.3],1,0,1);
    set(p,'LineStyle','--')
   export_fig([fig_path 'edgeInFar.pdf']) 
end








plot(radiusX,radiusY,'m','Linewidth',1.5);

% stinitko
bbx2 = [-7 -4 15 25];


plot(l1(1,:),l1(2,:),'k','Linewidth',1.5); plot(l2(1,:),l2(2,:),'k','Linewidth',1.5);



























