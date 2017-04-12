xlim = [-6 13];

t = xlim(1):0.05:xlim(2);
n_x = length(t);
sigma = 2;
mean1 = 0;
mean2 = 6.5; 

n_vrstev = 4;
vrstva = 3;

G1 = gaussMean(t,sigma,mean1)/1.5;
G2 = gaussMean(t,sigma,mean2);
G = G1+G2;
figure(1);clf
axis off


plot(t,G1,'g--'); hold on
plot(t,G2,'b--'); hold on
plot(t,G,'k','Linewidth',1.5); hold on
plot([t(1) t(end)],[0 0],'k-')


set(gca,'Xlim',xlim)

f_max = max(G);
set(gca,'Ylim',[0 f_max])

vrstvy = 0:f_max/(n_vrstev) :f_max;
vrstvy_y = zeros(n_vrstev,n_x);

idx_r = [];

vrstvy_y(vrstva-1,:) = vrstvy(vrstva)*ones(1,n_x);
greater = G >  vrstvy_y(vrstva-1,:);
x_fill = [];
y_fill = [];

for r = 1:n_x
    if greater(r)
        x_fill(end+1) = t(r);
        y_fill(end+1) = G(r);
        idx_r(end+1) = r;
    else
        if ~isempty(idx_r)

          	while ~isempty(idx_r)
               y_fill(end+1) = vrstvy_y(vrstva-1,idx_r(end));
               x_fill(end+1) = t(idx_r(end));
               idx_r(end) = [];
            end
%             fill(x_fill,y_fill,[0.7 0.7 0.7])
            x_fill = [];
            y_fill = [];
        end
    end
end

text_vrstvy = {'e','d','c','b','a'};
n_layer = 1;
for w = 1:(length(vrstvy)-1)
   vrstvy_y(w,:) = vrstvy(w)*ones(1,n_x);
   greater = G >  vrstvy_y(w,:);
   p1 = [];
   for r = 1:n_x
      if greater(r)
          if isempty(p1)
              p1 = t(r);
              text(p1-0.75,vrstvy(w),['\textbf{' text_vrstvy{n_layer} '}'],'Color','r','Interpreter','latex','FontSize',16)
              n_layer = 1+n_layer;
          else
              p2 = t(r);
              plot([p1 p2],[vrstvy(w) vrstvy(w)],'r','Linewidth',1.5);
              
              p1 = p2;
          end
      else
         p1 = []; 
      end
   end
   
end

% t1 = [mean1-1.5 0.123 1];
% t2 = [mean2+1.5+0.5 0.153 0];
% 
% q1 = [mean1 0.12 1];
% q2 = [mean2 0.15 0];
% 
% alpha = 0.8*0.5/norm(t1-q1);  % Size of arrow head relative to the length of the vector
% beta  = 0.007;  % Width of the base of the arrow head relative to the length
% vectarrow(t1,q1,'k', [1 1 0],alpha, beta)
% 
% text(t1(1)-2.5,t1(2),'1.peak','Interpreter','latex','fontsize',12)
% 
% 
% alpha = 0.8*0.5/norm(t2-q2);  % Size of arrow head relative to the length of the vector
% beta  = 0.007;  % Width of the base of the arrow head relative to the length
% vectarrow(t2,q2,'k', [1 1 0],alpha, beta)
% 
% text(t2(1)+0.1,t2(2),'2.peak','Interpreter','latex','fontsize',12)

axis off
set(gcf,'Position', [403 338 560 328]) 
set(gcf,'Color','w')


h = legend('$y_1 = k1\cdot e^{-\frac{(x-\mu_1)^2}{2\,\sigma}}$','$y_2 = k_2\cdot e^{-\frac{(x-\mu_2)^2}{2\,\sigma}}$','$y_1+y_2$');
set(h,'Interpreter', 'latex')
set(h,'Position', [0.1277 0.7559 0.2556 0.1821])
set(h,'Box','off')
set(h,'FontSize',12)

export_fig(['C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\' 'gaussIntersection.pdf'])


%%
figure(2);
clf;
text_vrstvy = {'e','d','c','b','a'};
axis([-0.5 0.5 -0.02 1])
text(0,0,['\textbf{' text_vrstvy{1} '}'],'Color','r','Interpreter','latex','FontSize',24)
hold on
l = 0.2;
fi = pi/6+pi/26;
t1 = [.02 .05 0];
q1 = [.02 .25 0];
alpha = 0.08*0.5/norm(t1-q1);  % Size of arrow head relative to the length of the vector
beta  = 0.4;  % Width of the base of the arrow head relative to the length
vectarrow(t1,q1,'k', [1 1 0],alpha, beta)

axis equal
text(-0.005,.29, ['\textbf{' text_vrstvy{2} '}'],'Color','r','Interpreter','latex','FontSize',24)
hold on

t1 = [.02 .35 0];
q1 = [.02-l*sin(fi) 0.35+l*cos(fi) 0];
alpha = 0.08*0.5/norm(t1-q1);  % Size of arrow head relative to the length of the vector
beta  = 0.4;  % Width of the base of the arrow head relative to the length
vectarrow(t1,q1,'k', [1 1 0],alpha, beta)


text(.02-l*sin(fi)-0.02,0.35+l*cos(fi)+0.04,['\textbf{' text_vrstvy{3} '}'],'Color','r','Interpreter','latex','FontSize',24)
hold on

t1 = [.02 .35 0];
q1 = [.02+l*sin(fi) 0.35+l*cos(fi) 0];
alpha = 0.08*0.5/norm(t1-q1);  % Size of arrow head relative to the length of the vector
beta  = 0.4;  % Width of the base of the arrow head relative to the length
vectarrow(t1,q1,'k', [1 1 0],alpha, beta)

text(.02+l*sin(fi)-0.02,0.35+l*cos(fi)+0.04,['\textbf{' text_vrstvy{4} '}'],'Color','r','Interpreter','latex','FontSize',24)

t1 = [.02+l*sin(fi) 0.35+l*cos(fi)+0.10 0];
q1 = [.02+l*sin(fi) 0.35+l*cos(fi)+0.30 0];
alpha = 0.08*0.5/norm(t1-q1);  % Size of arrow head relative to the length of the vector
beta  = 0.4;  % Width of the base of the arrow head relative to the length
vectarrow(t1,q1,'k', [1 1 0],alpha, beta)

text(.02+l*sin(fi)-0.02,0.35+l*cos(fi)+0.34, ['\textbf{' text_vrstvy{5} '}'],'Color','r','Interpreter','latex','FontSize',24)

axis off
set(gcf,'color','w')

export_fig(['C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\' 'treeGauss.pdf'])




