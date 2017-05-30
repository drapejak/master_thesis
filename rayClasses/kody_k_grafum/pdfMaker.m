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

% eval(['stone =' groupes{1} ';']);
% stone=addactray(stone, 'top');
% stone=solveraysen(stone,max(n_bound));
latexText = {};  
%% optimalizovane a kriterialni 
latexText(end+1,1) = {'\subsection*{Optimize and criterial}'};
for q = 1:length(groupes)
    structGroupes = [];
    newGroupes = load([groupFolder groupes{q}]);
    try
        structGroupes = load([savedGroupFolder groupes{q}]); % zatim nema vyznam, myslenka je takova, ze bychom nacetli predchozi vysledky a pridali pouze nove skupiny svazku
    catch
       warning(['New set of groupes for ' groupes{q} ' is determined.' ]) 
    end
    newGroupes = newGroupes.choosen_groupes;    
    
    
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
    
    groupes_sorted = struct('paths',GRall,'bounds_cnt',num2cell(val),'index',groupIndeces');
    save([groupFolder2 groupes{q} '_groupes'], 'groupes_sorted')
    
    eval(['stone =' groupes{1} ';']);
    
    stone=addactray(stone, 'top');
    stone=solveraysen(stone,max(n_bound+1));
        
    for w = 1:length(newGroupes)
        GR = GRall{w};      
        latexText(end+1,1) = {['\subsection*{Group ' groupIndeces{w} '}']};
        rayIdx = selectrays2(stone,{{1 GR{1,:}}});
        if ~isempty(rayIdx)
            figure(w);clf
            drawmodel(stone,rayIdx,'NoNewFigure');
            view([-9 60])
            set(gcf,'Color',[1 1 1])
            axis off
            zoom(gcf, 2.2)
            myprint([savedFIGFolder 'group' groupIndeces{w}],'-dpdf',10,10)
%             print([savedFIGFolder 'group' groupIndeces{w}],'-dpdf')
            delete(gcf)
            latexText(end+1: (end+nMiss),1) = {''};
            latexText(end+1: (end+nMiss),1) = {''};
            latexText(end+1,1) = {'\begin{figure}[h!]'};    %#ok                  
            latexText(end+1,1) = {'\centering'};            %#ok
            latexText(end+1,1) = {['\includegraphics[width = 8cm]{group' groupIndeces{w} '.pdf}']};     %#ok
            latexText(end+1,1) = {['\caption{Paths of rays in stone for group ' groupIndeces{w} '}']};  %#ok
            latexText(end+1,1) = {['\label{table:FigGroup' groupIndeces{w} '}']};   %#ok                
            latexText(end+1,1) = {'\end{figure}'};          %#ok
            
        end
        if size(GR,1)>6
            GR = GR(1:6,:);
            GR(end+1,:) = {'\dots'};
        end
        clear input;
        % Now use this table as input in our input struct:
        input.data = GR;
        % LaTex table caption:
        input.tableCaption = ['Paths of rays namely for group ' groupIndeces{w} ];
        input.tablePlacement ='h!';
        % LaTex table label:
        input.tableLabel = ['TableGroup' groupIndeces{w}];
        input.transposeTable = 1;
        latexText(end+1: (end+nMiss),1) = {''};
%         latexText(end+1) = {['\section{' 'Group ' groupIndeces{w} '}']}; %#ok 
        % Now call the function to generate LaTex code:
        latexTextNew = latexTable(input);
        latexText(end+1:(end+length(latexTextNew))) = latexTextNew;
        latexText(end+1,1) = {'\newpage'};          %#ok
    end
    
end

%% kandidatske
groupFolder         = [folder '\groupDistribution\resultsCand\'];
groupFolder2        = [folder '\groupDistribution\results_sortedCand\'];
latexText(end+1,1) = {'\subsection*{Candidates}'};
for q = 1:length(groupes)
    structGroupes = [];
    newGroupes = load([groupFolder groupes{q}]);
    try
        structGroupes = load([savedGroupFolder groupes{q}]); % zatim nema vyznam, myslenka je takova, ze bychom nacetli predchozi vysledky a pridali pouze nove skupiny svazku
    catch
       warning(['New set of groupes for ' groupes{q} ' is determined.' ]) 
    end
    newGroupes = newGroupes.choosen_groupes;    
    
    
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
    
    groupes_sorted = struct('paths',GRall,'bounds_cnt',num2cell(val),'index',groupIndeces');
    save([groupFolder2 groupes{q} '_groupes'], 'groupes_sorted')
    
    eval(['stone =' groupes{1} ';']);
    
    stone=addactray(stone, 'top');
    stone=solveraysen(stone,max(n_bound+1));
        
    for w = 1:length(newGroupes)
        GR = GRall{w};      
        latexText(end+1,1) = {['\subsection*{Group ' groupIndeces{w} '}']};
        rayIdx = selectrays2(stone,{{1 GR{1,:}}});
        if ~isempty(rayIdx)
            figure(w);clf
            drawmodel(stone,rayIdx,'NoNewFigure');
            view([-9 60])
            set(gcf,'Color',[1 1 1])
            axis off
            zoom(gcf, 2.2)
            myprint([savedFIGFolder 'group' groupIndeces{w}],'-dpdf',10,10)
%             print([savedFIGFolder 'group' groupIndeces{w}],'-dpdf')
            delete(gcf)
            latexText(end+1: (end+nMiss),1) = {''};
            latexText(end+1: (end+nMiss),1) = {''};
            latexText(end+1,1) = {'\begin{figure}[h!]'};    %#ok                  
            latexText(end+1,1) = {'\centering'};            %#ok
            latexText(end+1,1) = {['\includegraphics[width = 8cm]{group' groupIndeces{w} '.pdf}']};     %#ok
            latexText(end+1,1) = {['\caption{Paths of rays in stone for group ' groupIndeces{w} '}']};  %#ok
            latexText(end+1,1) = {['\label{table:FigGroup' groupIndeces{w} '}']};   %#ok                
            latexText(end+1,1) = {'\end{figure}'};          %#ok
            
        end
        if size(GR,1)>6
            GR = GR(1:6,:);
            GR(end+1,:) = {'\dots'};
        end
        clear input;
        % Now use this table as input in our input struct:
        input.data = GR;
        % LaTex table caption:
        input.tableCaption = ['Paths of rays namely for group ' groupIndeces{w} ];
        input.tablePlacement ='h!';
        % LaTex table label:
        input.tableLabel = ['TableGroup' groupIndeces{w}];
        input.transposeTable = 1;
        latexText(end+1: (end+nMiss),1) = {''};
%         latexText(end+1) = {['\section{' 'Group ' groupIndeces{w} '}']}; %#ok 
        % Now call the function to generate LaTex code:
        latexTextNew = latexTable(input);
        latexText(end+1:(end+length(latexTextNew))) = latexTextNew;
        latexText(end+1,1) = {'\newpage'};          %#ok
    end
    
end

% save LaTex code as file
fid=fopen([latexFolder 'text.tex'],'w');
[nrows,ncols] = size(latexText);
for row = 1:nrows
    fprintf(fid,'%s\n',latexText{row,:});
end
fclose(fid);
display(['... your LaTex code has been saved as ''text.tex'' in ' latexFolder ]);
    

    
    



