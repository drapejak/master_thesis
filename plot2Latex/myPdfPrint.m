function myPdfPrint(h,filemane,left,bottom,right,top)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','manual','PaperUnits','Inches',...
    'PaperPosition',[left bottom pos(3) pos(4)],...
    'PaperSize',[pos(3)+right, pos(4)+top])
print(h,filemane,'-dpdf','-r0')

