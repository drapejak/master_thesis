function [alpha, beta, R,T]...
                        =countboundkoef(u,n,nr,reflectionface)

%function [alpha, beta, totalreflection,Rp,Rl,Tp,Tl,Fp,Fl]...
%                        =countboundkoef(u,n,nr,reflectionface);
%
%Tato funkce spocita koeficienty odrazu a lomu paprsku.
%
%Vstupni parametry:
%   u-vektor urcujici smer sireni paprsku
%   n-normalovy vektor facety(roviny ve ktere lezi rozhrani), smerujici
%       do poloprostoru odkud prichazi paprsek.
%   nr-relativni idex lomu na prechodu (n2/n1, kde n1 je index prostredi ze ktereho paprsek prichazi).
%   reflectionface- 1-pokud se jedna o odraz na zrcadle, 0 pokuj dle
%       pouze o pruchod mezi prostrenimy ruzne opticke hustoty
%
%Vystup:
%   aplha-uhel dopadu (a odrazu)
%   beta-uhel lomu (NaN pokud doslo k totalnimu odrazu)
%   R-struktura obsahujici koeficienty zmeny odrazeneho paprsku
%   T-struktura obsahujici koeficienty zmeny lomeneho paprsku
%
%   format struktut R,T:
%       X.Ap-zmena amplitudy slozky svetla polarizovaneho v rovine kolme
%           k rovine dopadu.
%       X.Al-zmena amplitudy slozky svetla polarizovaneho v rovine rovnobezne
%           s rovinou dopadu.
%       X.Fp-zmena faze svetla polarizovaneho v rovine kolme
%          k rovine dopadu pri odrazu. (v radianech)
%       X.Fl-zmena faze svetla polarizovaneho v rovine rovnobezne
%       s rovinou dopadu pri odrazu. (v radianech)
%

constants;

%alpha- uhel dopadu
alpha=acos((-u'*n)/(norm(u)*norm(n)));

if reflectionface
    %PRECHOD JE ODRAZ NA ZRCADLE
    beta=NaN;   %uhel lomu
    R=struct('Ap',-1,'Al',-1,'Fp',0,'Fl',0);
    %????????????????????????????????????????????????
    %nema dojit ke zmene faze ????? jako pri normalnik totalnim odrazu?
    %nebo jako prri prechodu do prostredi s n->inf ???
    %Fl=pi;%********
    %Fp=pi;%????????
    %????????????????????????????????????????????????
    % TODO opravit!!! Smutny
    T=NaN;
else
    %PRECHOD JE NORMALNI
    beta_kf=sin(alpha)/nr;

    if beta_kf>=1
        %DOJDE K TOTALNIMU ODRAZU
        beta=NaN;

        %zmena inensity a faze
        f_tmp= (sin(alpha)^2-nr^2)^(1/2) / cos(alpha);
        R=struct('Ap',1,'Al',1,'Fl',-2*atan(f_tmp/nr^2),'Fp',-2*atan(f_tmp));
        T=NaN;
    else
        %DOJDE K ROZDELENI PAPRSKU
        beta=asin(beta_kf);

        if (abs(alpha)<MIN_ANGLE),
            %paprsek dopada kolmo %????? zkontrolovat
            r_tmp=(nr-1)/(nr+1);
            R=struct('Ap',r_tmp,'Al',-r_tmp,'Fl',0,'Fp',0);
            t_tmp=2/(nr+1);
            T=struct('Ap',t_tmp,'Al',t_tmp,'Fl',0,'Fp',0);
        else
            %paprsek nedopada kolmo
            R=struct(...
                'Ap',-sin(alpha-beta)/sin(alpha+beta), ...
                'Al',tan(alpha-beta)/tan(alpha+beta), ...
                'Fl',0,'Fp',0 );
            T=struct(...
...%                'Ap',2*sin(beta)*cos(alpha)/sin(alpha+beta), ...
...%                'Al',2*sin(beta)*cos(alpha)/(sin(alpha+beta)*cos(alpha-beta)), ...
...% Opraveno? Smutny 2015/3/17, vysledek se lisil od Golsteina a amplituda prvku muellerovy matice byla vetsi nez 1
                'Ap',sqrt(cos(beta)/cos(alpha))*2*sin(beta)*cos(alpha)/sin(alpha+beta), ...
                'Al',sqrt(cos(beta)/cos(alpha))*2*sin(beta)*cos(alpha)/(sin(alpha+beta)*cos(alpha-beta)), ...
                'Fp',0,'Fl',0);
            %pozn. neni treba osetrovat pripad alpha+beta=pi/2  matlab umi poctitat s nekonecnem
         end;
    end;
end;


% %%%TMP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isnan(beta) return; end;
%
% k=cos(beta)/cos(alpha)*nr;
%
% ETl= T.Al^2*k;
% ERl= R.Al^2;
% El=ETl+ERl;
%
% ETp= T.Ap^2*k;
% ERp= R.Ap^2;
% Ep=ETp+ERp;
%
% disp(sprintf('L: C[%f] R[%f] T[%f]',El,ERl,ETl));
% disp(sprintf('P: C[%f] R[%f] T[%f]',Ep,ERp,ETp));



