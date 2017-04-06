tvar = test_2_blunt_13_7_16(15);
pathProtocol = 'C:\Users\kubin\Desktop\GitFolder\doc_thesis\figures';
drawmodel(tvar)
axis off
set(gcf,'Color','w')
set(gcf,'Position', [680 672 587 306])
set(gca,'CameraPosition', [-9.5238 5.7718 2.2896])
set(gca,'CameraTarget', [0.6655 0.8487 0.9426])
set(gca,'CameraViewAngle', 1.1512)

export_fig([ pathProtocol filesep 'edge_ex.pdf'], '-q101')