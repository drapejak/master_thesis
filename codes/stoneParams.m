stone = viva12;
textSize = 12;
drawmodel(stone,[],'EdgeColor','k')
axis off
set(gcf,'Color','w')
facenames = getfacenames(stone);
for r = 1:12
    pnt = getfacecentroid(stone,r);
   text(pnt(1)-0.05,pnt(2),pnt(3),['\textbf{' facenames{r} '}'],'FontSize',textSize,'Interpreter','Latex') 
end

pnt = getfacecentroid(stone,13);
text(pnt(1)-0.05,pnt(2),pnt(3),['\textbf{' 'TOP' '}'],'FontSize',textSize,'Interpreter','Latex') 
export_fig (['C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\' 'vivi12Facets1.pdf'])

%%
stone =viva12model([], 1.514, 4,1.6, 1.6);
drawmodel(stone,[],'EdgeColor','k')
axis off
set(gcf,'Color','w')
view(90,0)
facenames = getfacenames(stone);
for r = [ 1 2 11 12]
    pnt = getfacecentroid(stone,r);
   text(pnt(1),pnt(2)-0.05,pnt(3),['\textbf{' facenames{r} '}'],'FontSize',textSize,'Interpreter','Latex') 
end

for r = [ 10]
    pnt = getfacecentroid(stone,r);
%    text(pnt(1),pnt(2)-0.09,pnt(3),facenames{r},'FontSize',12,'Interpreter','Latex') 
   text(pnt(1),pnt(2)-0.24,pnt(3),'$\mathbf{h}$' ,'FontSize',textSize,'Interpreter','Latex') 
end

% for r = [3]
%     pnt = getfacecentroid(stone,r);
%    text(pnt(1),pnt(2)-0.02,pnt(3),facenames{r},'FontSize',12,'Interpreter','Latex') 
% end
r = 13; % top
pnt = getfacecentroid(stone,r);
text(pnt(1),pnt(2)-0.05,pnt(3)+0.03,['\textbf{' 'TOP' '}'],'FontSize',textSize,'Interpreter','Latex') 
text(pnt(1),pnt(2)-0.05,pnt(3)+0.12,'$\mathbf{d_{TOP}}$','FontSize',textSize,'Interpreter','Latex') 

r = 14; % bottom
pnt = getfacecentroid(stone,r);
text(pnt(1),pnt(2)-0.05,pnt(3)-0.04,['\textbf{' 'BOT' '}'],'FontSize',textSize,'Interpreter','Latex') 
text(pnt(1),pnt(2)-0.05,pnt(3)-0.12,'$\mathbf{d_{BOT}}$','FontSize',textSize,'Interpreter','Latex') 

% kota h
plot3([0 0],[-0.205 -0.575],[0 0],'k','Linewidth',1.5)
plot3([0 0],[-0.5 -0.575],[-0.375 -0.375],'k','Linewidth',1.5)

p1 = [0 -0.55 0];
p2 = [0 -0.55 -0.25];
p3 = [0 -0.55 -0.375];

% myQuiver(p2,p1,0.4)
% myQuiver(p2,p3,0.8)
vectarrow(p2,p1,'k',[0 1 1])
vectarrow(p2,p3,'k',[0 1 1])


sy = -0.375;
lkota = 0.110;
dkota = 0.025;

% kota d Bottom
plot3([0 0],[0.5 0.5],[sy   sy-lkota],'k','Linewidth',1.5)
plot3([0 0],[-0.5 -0.5],[sy sy-lkota],'k','Linewidth',1.5)

p1 = [0 -0.5 sy-lkota+dkota];
p2 = [0 0    sy-lkota+dkota];
p3 = [0 0.5  sy-lkota+dkota];

% q1 = myQuiver(p2,p1,0.4);
% q2 = myQuiver(p2,p3,0.8);

vectarrow(p2,p1,'k',[0 1 1])
vectarrow(p2,p3,'k',[0 1 1])


sy = 0;
lkota = 0.110;
dkota = 0.025;
dy = 0.205;
% kota d top
plot3([0 0],[dy dy],[sy   sy+lkota],'k','Linewidth',1.5)
plot3([0 0],[-dy -dy],[sy sy+lkota],'k','Linewidth',1.5)

p1 = [0 -dy sy+lkota-dkota];
p2 = [0 0   sy+lkota-dkota];
p3 = [0 dy  sy+lkota-dkota];

% q1 = myQuiver(p2,p1,0.4);
% q2 = myQuiver(p2,p3,0.8);

vectarrow(p2,p1,'k',[0 1 1])
vectarrow(p2,p3,'k',[0 1 1])

% kota lem
plot3([0 0],-[-0.5 -0.575],[-0.3 -0.3],'k','Linewidth',1.5)
plot3([0 0],-[-0.5 -0.575],[-0.375 -0.375],'k','Linewidth',1.5)

p1 = [0 0.55 -0.3];
p2 = [0 0.55 -0.225];
p3 = [0 0.55 -0.375];
p4 = [0 0.55 -0.45];
% myQuiver(p2,p1,0.4)
% myQuiver(p2,p3,0.8)
vectarrow(p2,p1,'k',[0 1 1])
vectarrow(p4,p3,'k',[0 1 1])
plot3([p1(1) p3(1)],[p1(2) p3(2)],[p1(3) p3(3)],'k','Linewidth',1.5)

text(p2(1),p2(2)+0.02,p2(3)-0.05,'$\mathbf{h_{RF}}$','FontSize',textSize,'Interpreter','Latex') 

export_fig (['C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\' 'vivi12Facets2.pdf'])

%%






