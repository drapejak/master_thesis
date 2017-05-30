autoOptimization;
child  = get(gcf,'Children');

for w = 1:8
    figure
    set(child(w),'Parent',gcf);
    set(child(w),'Units','pixels')
    pos = get(child(w),'Position');
    pos2 = get(gcf,'Position');
    set(gcf,'InvertHardcopy','off')
    set(child(w),'Visible', 'on');
     set(child(w),'FontSize', 8);
    set(gcf,'Position',[pos2(1:2) pos(3:4) ])
    set(child(w),'Position',[0 0 pos(3:4) ])
    j = 1;   
    for i = 1: 80000000
        j = j+1;
   end
    printType = '-dpdf';
    changeAxes = 0;

    sizx = pos(3)/100;     sizy = pos(4)/100;
    posoffx = 0;    posoffy = 0;
    zoomi = 1;
    myprint(sprintf('%d',w),printType,sizx,sizy,zoomi,posoffx,posoffy,changeAxes)
    export_fig(sprintf('%d2.pdf',w),gcf)
end