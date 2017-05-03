%% 19/4/2017 odraz sv?tla v klínu
meanThres = -1;  % ocekavame, ze pixely, ktere jsou mensi, nez meanThres*mean(im) jsou mimo kouli
backgroundParam = 201;  % velikost filtru a zaroven velikost sigma pro odhad pozadi

backgroundk = 0.3;

folder= 'D:\CuStoRe\20160914-viva\';
image2 = imread([folder 'viva_15_3' '.tiff']);
rays_detected2 = detectStoneFaster(image2,0);

% imf = convolve2(double(image2),fspecial('gaussian',filterSize,filterSigma) ,'same');
im_mean = meanThres*mean(image2(:));

% Odfiltrovani pozadi od obrazku
im_background = convolve2(image2,fspecial('gaussian',backgroundParam,backgroundParam) ,'wrap')+ backgroundk*im_mean;
Image_filtered = max(double(image2) - im_background,0);

punch2 = [763 1370 780 1330];
figure(233);
plotImage(Image_filtered,1,0); hold on
axis(punch2)
axis off
set(gcf,'Color','w')
set(gcf,'Position',[262 26 894 657])


mask = (rays_detected2.Position(:,1) > punch2(1)) &...
    (rays_detected2.Position(:,1) < punch2(2))& ...
    (rays_detected2.Position(:,2) > punch2(3))& ...
    (rays_detected2.Position(:,2) < punch2(4));

rays_choosen = find(mask);
px_width = 12;
vyber = [53 51 50 48 72 40 18 7 5 12 28 85 92 89 56 79  71 37 16 6 4 9 27 83 91 88 55 78];
names = ['1B' '3B' '5E' '7H' repmat({'3A'},1,12) repmat({'5D'},1,12)];
my_color = ['r','r','r','r',repmat('b',1,12),repmat('m',1,12)];
for r  = 1:length(vyber)
   pnt =  rays_detected2.Position(rays_choosen(vyber(r)),:);
   plot(pnt(1),pnt(2),[my_color(r) 'o'],'Linewidth',1,'MarkerSize',12)
   text(pnt(1)+12,pnt(2),names{r},'FontSize',18,'Color',my_color(r))


end
%   export_fig(['C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures\' 'wedge_example' '.pdf'])
