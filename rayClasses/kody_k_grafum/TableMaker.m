% prohleda urcene skupiny svazku, porovna s ulozenymi skupinami svazku a
% nove skupiny zanese do latexovskeho souboru i cestou 
% spustit a vybrat Add t
folder = fileparts(which('campaths.m'));

groupFolder         = [folder '\groupDistribution\results\'];
groupFolder2        = [folder '\groupDistribution\results_sorted\'];
latexFolder         = [folder '\Dokumentace\SkupinySvazku\'];
savedGroupFolder    = [latexFolder 'groupes\'];
savedFIGFolder      = [latexFolder 'figures\'];

addpath([latexFolder 'latexTable'])

groupes  = getAllFiles(groupFolder,'.mat');
acbde = {'A','B','C','D','E','F','G','H','I','J','K','L',...
        'M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

% eval(['stone =' groupes{1} ';']);
% stone=addactray(stone, 'top');
% stone=solveraysen(stone,max(n_bound));
latexText = {};  

%% hodnoty do tabulky
% optimalizovane a kriterialni 
groupFolder         = [folder '\groupDistribution\results\'];
abcde_idx = ones(1,11);
structGroupes = [];
q = 1;
    newGroupes = load([groupFolder groupes{q}]);
    try
        structGroupes = load([savedGroupFolder groupes{q}]); % zatim nema vyznam, myslenka je takova, ze bychom nacetli predchozi vysledky a pridali pouze nove skupiny svazku
    catch
       warning(['New set of groupes for ' groupes{q} ' is determined.' ]) 
    end
    newGroupes = newGroupes.choosen_groupes;
    % pridani nove skupiny svazku
    newGroupes(end+1) = newGroupes(1);
    for w = 1:12
       newGroupes{end}{w}([2,4]) = newGroupes{end}{w}([4,2]);
    end
    
    % nacteme si vysledky a seradime podle odrazu 
    nMiss = 3;
    GRall = {};
    for w = 1:length(newGroupes)
        GR = {};
        for r = 1:length(newGroupes{w})
            GR(r,:) = newGroupes{w}{r}(2:end);
        end
        GRall{w} = GR;
    end    
    
    % kandidatske 
    groupFolder         = [folder '\groupDistribution\resultsCand\'];
    structGroupes = [];
    newGroupes = load([groupFolder groupes{q}]);
    try
        structGroupes = load([savedGroupFolder groupes{q}]); % zatim nema vyznam, myslenka je takova, ze bychom nacetli predchozi vysledky a pridali pouze nove skupiny svazku
    catch
       warning(['New set of groupes for ' groupes{q} ' is determined.' ]) 
    end
    newGroupes = newGroupes.choosen_groupes;    
    
    
    % nacteme si vysledky a seradime podle odrazu 
    
    for w = 1:length(newGroupes)
        GR = {};
        for r = 1:length(newGroupes{w})
            GR(r,:) = newGroupes{w}{r}(2:end);
        end
        GRall{end+1} = GR;
    end   
    
    
    [~, ncols]  = cellfun(@size,GRall);
    [val,idx]   = sort(ncols);
    n_cnt   = accumarray(val',ones(length(idx),1));
    n_cnt   = n_cnt(n_cnt>0)';
    n_bound = unique(val);
    GRall = GRall(idx);
    groupIndeces = cell(sum(n_cnt),1);
    kk = 1;
    for k = 1:length(n_bound)
       for l = 1: n_cnt(k)
           groupIndeces(kk) = {[int2str(n_bound(k)) int2str(l-1)]};
           kk = kk+1;
       end
    end
    class_index = zeros(1,sum(n_cnt));
    for k = 1:sum(n_cnt)
        class_index(k) =  abcde_idx(val(k));       
        abcde_idx(val(k)) = abcde_idx(val(k)) +1;
    end
    
    % prerovnani podle 1. fasety
    for w = 1:length(newGroupes)
        GR = GRall{w};
        gr_length = cellfun(@length,GR(:,1));
        mask = find(gr_length == min(gr_length));
        firstFacet = mask(~cellfun(@isempty, regexp(GR(mask,1),'1')));
        if ~isempty(firstFacet)
            GR = circshift(GR, [-firstFacet+1, 0]);
        end
        GRall{w} = GR;
    end
    [len,~] = cellfun(@size,GRall);
    
    groupes_sorted = struct('paths',GRall,'bounds_cnt',num2cell(val),'index',groupIndeces','len',num2cell(len),'class_index',num2cell(class_index));
    save([groupFolder2 groupes{q} '_groupesAll'], 'groupes_sorted')

    
%% tvorba obrazku 
    
%     eval(['stone =' groupes{1} ';']);
%     
%     stone=addactray(stone, 'top');
%     stone=solveraysen(stone,max(n_bound+1));
latexText = {};
viva_classes = [];
S = load('D:\CuStoRe\20160914-viva\20162209_resultsAuto\viva_11_2Optim.mat');
close all;
try
    stone = S.dataOptim.model.weray;
catch
    stone = S.dataOptim.model;
end
    for w = 1:length(GRall)
        GR = GRall{w};      
        latexText(end+1,1) = {['\subsection*{Group ' groupIndeces{w} '}']};
        en_max  = 0;
        rayIdx = [];
        for ww = 1:groupes_sorted(w).len    
           rayIdxAct = selectrays2(stone,{{1 GR{ww,:}}});
           GR(ww) = {{1 GR{ww,:}}};
           if ~isempty(rayIdxAct)
               energyAct = getrayenergy(stone,rayIdxAct);
               if energyAct > en_max
                    en_max =   energyAct;
                    rayIdx = rayIdxAct;
               end
           end
        end
        className = [num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index}];
        if strcmp(className,'3B')
            className = '3C';
        elseif strcmp(className,'3C')
            className = '3B';
        end
        eval(['viva_classes.c_' className ' = GR(:,1);'])
        
        
%         if ~isempty(rayIdx)
%             figure(w);clf
%             drawmodel(stone,rayIdx,'NoNewFigure','ExpP','AxisUp',0.1);
%             view([-9 60]);
%             set(gcf,'Color',[1 1 1]);
%             axis off;
%             zoom(gcf, 2.2);
% %             export_fig([savedFIGFolder 'group' num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index} '.pdf' ])
%             myprint([savedFIGFolder 'group' num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index} ],'-dpdf',10,10)
% %             print([savedFIGFolder 'group' groupIndeces{w}],'-dpdf')
%             
%         else 
%             warning(['Missing ' num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index}])
%         end
    end
 save('viva_classes','viva_classes')          

%%
latexText = {};  
n_paths = 8;
for w = 1:length(groupes_sorted)
   for r = 1:groupes_sorted(w).bounds_cnt
      class = '';
      rcount = '';
      if r == 1
          pathText = '';
           class =  ['\textbf{' num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index} '}'];
           rcount = num2str(groupes_sorted(w).len);
      end
      
      
      if groupes_sorted(w).len > 1
          pathText = ' & ';
          for rr = 1:n_paths
             pText = groupes_sorted(w).paths{rr,r};
             if strcmp(pText,'bottom')
                 pText = 'BOT';
             elseif strcmp(pText,'top')
                 pText = 'TOP';
             end
             pathText = [pathText pText ' & '];
          end
          pathText=  [pathText '$\dots$'];
          latexText(end+1,1) = {[class pathText ' & ' rcount '\\']};
      else
          pText = groupes_sorted(w).paths{1,r};
             if strcmp(pText,'bottom')
                 pText = 'BOT';
             elseif strcmp(pText,'top')
                 pText = 'TOP';
             end
          if r == groupes_sorted(w).bounds_cnt
              class =  ['\textbf{' num2str(groupes_sorted(w).bounds_cnt) acbde{groupes_sorted(w).class_index} '}'];
              rcount = num2str(groupes_sorted(w).len);
              pathText = [' & \multicolumn{' int2str(n_paths+1) '}{l}{' pathText   pText '} \vline '];
              latexText(end+1,1) = {[class pathText ' & ' rcount '\\']};        
          else
              pathText = [pathText   pText '-'];
          end

          
      end
      
      
      
      
   end
   latexText(end+1,1) = {'\hline \hline'};
    
end
%%
% save LaTex code as file
fid=fopen([latexFolder 'text.tex'],'w');
[nrows,ncols] = size(latexText);
for row = 1:nrows
    fprintf(fid,'%s\n',latexText{row,:});
end
fclose(fid);
display(['... your LaTex code has been saved as ''text.tex'' in ' latexFolder ]);


    



