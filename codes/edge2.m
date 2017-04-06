addpath('D:\TDV\tdv_toolbox')
addpath C:\Users\kubin\Desktop\GitFolder\doc_thesis\codes
fig_path = 'C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\';

% definice bbx
bbx = [5 -4 15 4];

p1  = [bbx(1);bbx(4);1];
p2  = [bbx(3);bbx(4);1];
p3  = [bbx(1);bbx(2);1];
p4  = [bbx(3);bbx(2);1];
bbx_points = [p1 p2 p4 p3 p1]; 
figure
% plot(bbx_points(1,:),bbx_points(2,:),'b')
 hold on
 axis equal

% bbx lines 
top_line= cross(p1,p2);
bottom_line = cross(p3,p4);
l_line  = cross(p1,p3);
r_line  = cross(p2,p4);

%----------------------------%
%----------------------------%
%----------------------------%
addpath('D:\TDV\tdv_toolbox')

fi = 7/10*pi; % pi/2 - cca pi
radius_l = 3;
samples = 150;
vector = [0 1 0]'; % rovna faseta


% volba indexu lomu zmeni pripad lomu na hrane 
no = 1.0;
ni = 1.5;
fromGlass = ni > no;

if fromGlass
    vi = [0 1 0]';
else
    vi = [1 0.1 0]';
    vi = [1 0.05 0]';
    vi = [0.1 1 0]';
    vi = [0 -1 0]';
end



% prunik faset
xmax = 11;
ymin = 1;


ref_fi = pi - fi;
intersPoint = [xmax ymin]';

l1 = [bbx(1) xmax-radius_l;
      ymin  ymin];
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

plot(radiusX,radiusY,'k','Linewidth',1.5)


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
            l1_points = [cross(rayLine(:,w),top_line) cross(rayLine(:,w),r_line)];
        end
    l1_pointsE = p2e(l1_points);
    out =  logical(isinbbx(l1_pointsE, bbx));
    rayPoint(:,w) = l1_points(:,find(out,1,'first'));
    rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
%     plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-b')
    myq =  quiver(rayPoint(1,w),rayPoint(2,w),xEdge(w)-rayPoint(1,w),yEdge(w)-rayPoint(2,w),'b','AutoScale','off');
    if length(legsource) < 2
        legsource = [legsource myq];
    end
    else
        maskOk(w) = false;
    end
end
mask = mask(maskOk);

rayLine   = zeros(3,length(mask));
rayPoint = zeros(3,length(mask));
fio  = zeros(1,length(mask)); 
fio2 = zeros(1,length(mask)); 

cosAlfa = zeros(1,length(mask)); % uhel odrazu
cosBeta = zeros(1,length(mask)); % uhel lomu
for w = mask
    vo = real(ray_reflection(n_f(:,w),vi));
    vectOut = edge_points(:,w) - 20*vo;
    fio(w) = atan2(vo(1),vo(2));
    cosAlfa(w)   =  vo'*n_f(:,w)/(norm(vo)*norm(n_f(:,w)));
    
    rayLine(:,w) = cross(vectOut,edge_points(:,w));
    rayLine(:,w) = rayLine(:,w)/rayLine(3,w);
    
    if fromGlass
        l1_points = [cross(rayLine(:,w),bottom_line) cross(rayLine(:,w),l_line)];
    else
        l1_points = [cross(rayLine(:,w),top_line) cross(rayLine(:,w),r_line)];
    end
    l1_pointsE = p2e(l1_points);
    out =  logical(isinbbx(l1_pointsE, bbx));
    rayPoint(:,w) = l1_points(:,find(out,1,'last'));
    rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
%     plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-r')
    myq =  quiver(xEdge(w),yEdge(w),-xEdge(w)+rayPoint(1,w),-yEdge(w)+rayPoint(2,w),'r','AutoScale','off');
    
    if length(legsource) < 3
        legsource = [legsource myq];
    end
    
    vo = real(ray_refraction(n_f(:,w),vi,ni,no));
    cosBeta(w)   =  vo'*-n_f(:,w)/(norm(vo)*norm(-n_f(:,w)));
    
        
        fio2(w) = atan2(vo(1),vo(2));

        vectOut = edge_points(:,w) + 20*vo;

        rayLine(:,w) = cross(vectOut,edge_points(:,w));
        rayLine(:,w) = rayLine(:,w)/rayLine(3,w);

        if fromGlass
           l1_points = [cross(rayLine(:,w),top_line) cross(rayLine(:,w),r_line) cross(rayLine(:,w),l_line)];
        else
           l1_points = [cross(rayLine(:,w),bottom_line) cross(rayLine(:,w),r_line)];
        end
        l1_pointsE = p2e(l1_points);
        out =  logical(isinbbx(l1_pointsE, bbx));
        rayPoint(:,w) = l1_points(:,find(out,1,'first'));
        rayPoint(:,w) = rayPoint(:,w)/rayPoint(3,w);
    %     plot([rayPoint(1,w),xEdge(w)],[rayPoint(2,w),yEdge(w)],'-r')
    if abs(cosBeta(w))-1e-10 > 0 % totalni obraz
        myq = quiver(xEdge(w),yEdge(w),-xEdge(w)+rayPoint(1,w),-yEdge(w)+rayPoint(2,w),'g','AutoScale','off');
        if length(legsource) < 4
            legsource = [legsource myq];
        end
    end
end

h = legend(legsource,'edge','input light','reflected','refracted');
set(h,'Interpreter','latex')
if samples < 50
    if fromGlass
    %     text(13,-4,'$n_1 = 1.5$','Interpreter','Latex')
    text(13,-3.5,'Glass','Interpreter','Latex','FontSize',12)
        text(13, 0,'Air','Interpreter','Latex','FontSize',12) 
        
        export_fig([fig_path 'edgeOut.pdf'])
    else
       text(5.5,0,'Glass','Interpreter','Latex','FontSize',12)
        text(5.5, 2,'Air','Interpreter','Latex','FontSize',12)
       set(h,'Position', [0.1500 0.2104 0.1994 0.1546])
       export_fig([fig_path 'edgeIn.pdf']) 
    end
end
% fresnelovy rovnice odrazu
Rs  = (ni*cosAlfa - no*cosBeta)./(ni*cosAlfa + no*cosBeta);
Rp  = (ni*cosBeta - no*cosAlfa)./(ni*cosBeta + no*cosAlfa);
R = (Rs.^2+Rp.^2)/2;
T = 1-R;
hustota = [1 min(abs(diff(fio2)))./abs(diff(fio2))];

S_plocha = zeros(1,length(mask));
for w = mask
   S_plocha(w) = vi'*n_f(:,w)/(norm(vi)*norm(n_f(:,w)));   
end

if samples > 100
    figure;
    plot(fio,abs(S_plocha).*R,'r','Linewidth',1.5)
    set(gca,'Xlim',[min(fio) max(fio)])
    set(gca,'Ylim',[0 1])
    grid on
    ylabel('Reflection coefficient', 'Interpreter', 'latex')
    xlabel('$\varepsilon \left[ rad \right]$', 'Interpreter', 'latex')
    set(gca,'XTickMode','manual')
    set(gcf,'Position', [680 702 378 276])
    set(gcf,'Color',[1 1 1])
    if fromGlass
        set(gca,'XTick', [-pi -7/8*pi -3/4*pi -5/8*pi -pi/2 -3/8*pi -pi/4 ])

        set(gca,'TickLabelInterpreter','latex')
        set(gca,'XTickLabel',{'$-\pi$','$-\frac{7}{8}\pi$','$-\frac{3}{4}\pi$',...
            '$-\frac{5}{8}\pi$','$-\frac{\pi}{2}$','$-\frac{3}{8}\pi$','$-\frac{\pi}{4}$'})
        %'edgeIn_reflection'.eps rucne
            export_fig([fig_path 'edgeOut_reflection.pdf'])
    else
        set(gca,'XTick',  0:pi/8:3/4*pi)

        set(gca,'TickLabelInterpreter','latex')
        set(gca,'XTickLabel',{'$0$','$\frac{\pi}{8}$','$\frac{\pi}{4}$',...
            '$\frac{3}{8}\pi$','$\frac{\pi}{2}$','$\frac{5}{8}\pi$','$\frac{3}{4}\pi$'}) 
        %'edgeIn_reflection'.eps rucne
            export_fig([fig_path 'edgeIn_reflection.pdf'])
    end
    

    refrMask = T>1e-10;
    figure;
    plot(fio2(refrMask),abs(S_plocha(refrMask)).*T(refrMask).*hustota(refrMask),'g','Linewidth',1.5)
    set(gca,'Xlim',[min(fio2(refrMask)) max(fio2(refrMask))])
    set(gca,'Ylim',[0 1])
    grid on
    ylabel('Refraction coefficient', 'Interpreter', 'latex')
    xlabel('$\varepsilon \left[ rad \right]$', 'Interpreter', 'latex')
    set(gca,'XTickMode','manual')
    set(gcf,'Position', [680 702 378 276])
    set(gcf,'Color',[1 1 1])
    if fromGlass
        set(gca,'XDir','reverse')
        set(gca,'XTick', [-pi/4 -3/16*pi -2/16*pi -1/16*pi ])

        set(gca,'TickLabelInterpreter','latex')
        set(gca,'XTickLabel',{'$-\frac{\pi}{4}$','$-\frac{3}{16}\pi$','$-\frac{1}{8}\pi$',...
        '$-\frac{1}{16}\pi$'})
        %'edgeOut_refraction'.eps rucne
        export_fig([fig_path 'edgeOut_refraction.pdf'])
    else
        set(gca,'XTick', [-pi -31/32*pi -15/16*pi -29/32*pi])

        set(gca,'TickLabelInterpreter','latex')
        set(gca,'XTickLabel',{'$-\pi$','$-\frac{31}{32}\pi$','$-\frac{15}{16}\pi$','$-\frac{29}{32}\pi$'}) 
        %'edgeOut_refraction'.eps rucne
        export_fig([fig_path 'edgeIn_refraction.pdf'])
    
    end
    
    
end




























