% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function CoincidenceAngle_Anisotropic(Hybrid,Crop,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,PNGresolution,SColor,AColor,BColor,A,ALamb,AntisymmetricModes,AShear,AScholte,B,BLamb,BoxLineWidth,BShear,BScholte,Couplant,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeModeLabels,HeadLine,HigherOrderModes,LambModes,LineWidth,Material,ModeLabels,ModeLabel1X,ModeLabel2X,ModeLabel3X,PDF,FileName,PlateThickness,PNG,PropagationAngle,S,ShearHorizontalModes,ScholteModes,SLamb,SShear,SScholte,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled)
%#ok<*FXUP>
%#ok<*AGROW>
%#ok<*CHAIN>
f = figure('Name','Coincidence angle','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
hold on
if  Symmetric
    if  ~Decoupled
        if  SymmetricModes && ~isempty(S{1})
            S{1}(:,5) = real(asind(Couplant.Velocity/1e3./S{1}(:,4))); % Snell's law of refraction
            s = plot(S{1}(:,XAxisMode),S{1}(:,5),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(S{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(S{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)),S{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,S{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,S{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            S{2}(:,5) = real(asind(Couplant.Velocity/1e3./S{2}(:,4)));
            s(2) = plot(S{2}(:,XAxisMode),S{2}(:,5),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(S{2}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(S{2}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)),S{2}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,S{2}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,S{2}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 3:size(S,2)
                    S{i}(:,5) = real(asind(Couplant.Velocity/1e3./S{i}(:,4)));
                    s(i) = plot(S{i}(:,XAxisMode),S{i}(:,5),'LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  AntisymmetricModes && ~isempty(A{1})
            A{1}(:,5) = real(asind(Couplant.Velocity/1e3./A{1}(:,4)));
            a = plot(A{1}(:,XAxisMode),A{1}(:,5),'LineWidth',LineWidth,'Color',AColor);
            if  ModeLabels
                z = find(abs(A{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(A{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)),A{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,A{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,A{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(A,2)
                    A{i}(:,5) = real(asind(Couplant.Velocity/1e3./A{i}(:,4)));
                    a(i) = plot(A{i}(:,XAxisMode),A{i}(:,5),'LineWidth',LineWidth,'Color',AColor);
                end                
            end                
        end
    else
        if  LambModes && SymmetricModes && ~isempty(SLamb{1})
            SLamb{1}(:,5) = real(asind(Couplant.Velocity/1e3./SLamb{1}(:,4)));
            sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,5),'LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(SLamb{1}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(SLamb{1}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)),SLamb{1}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,SLamb{1}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,SLamb{1}(z(1),5)-(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(SLamb,2)
                    SLamb{i}(:,5) = real(asind(Couplant.Velocity/1e3./SLamb{i}(:,4)));
                    sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,5),'LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
            ALamb{1}(:,5) = real(asind(Couplant.Velocity/1e3./ALamb{1}(:,4)));
            aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,5),'LineWidth',LineWidth,'Color',AColor);
            if  ModeLabels
                z = find(abs(ALamb{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(ALamb{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)),ALamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,ALamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,ALamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(ALamb,2)
                    ALamb{i}(:,5) = real(asind(Couplant.Velocity/1e3./ALamb{i}(:,4)));
                    aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,5),'LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
            SShear{1}(:,5) = real(asind(Couplant.Velocity/1e3./SShear{1}(:,4)));
            sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            if  ModeLabels
                z = find(abs(SShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(SShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)),SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 2:size(SShear,2)
                    SShear{i}(:,5) = real(asind(Couplant.Velocity/1e3./SShear{i}(:,4)));
                    sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
            for i = 1:size(AShear,2)
                AShear{i}(:,5) = real(asind(Couplant.Velocity/1e3./AShear{i}(:,4)));
                aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        SScholte{1}(:,5) = real(asind(Couplant.Velocity/1e3./SScholte{1}(:,4)));
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                SScholte{i}(:,5) = real(asind(Couplant.Velocity/1e3./SScholte{i}(:,4)));
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        AScholte{1}(:,5) = real(asind(Couplant.Velocity/1e3./AScholte{1}(:,4)));
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                AScholte{i}(:,5) = real(asind(Couplant.Velocity/1e3./AScholte{i}(:,4)));
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end    
else
    if  ~Decoupled
        if  ~isempty(B{1})
            B{1}(:,5) = real(asind(Couplant.Velocity/1e3./B{1}(:,4)));
            b = plot(B{1}(:,XAxisMode),B{1}(:,5),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(B{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)),B{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,B{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,B{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            B{2}(:,5) = real(asind(Couplant.Velocity/1e3./B{2}(:,4)));
            b(2) = plot(B{2}(:,XAxisMode),B{2}(:,5),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{2}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(B{2}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)),B{2}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,B{2}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,B{2}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            B{3}(:,5) = real(asind(Couplant.Velocity/1e3./B{3}(:,4)));
            b(3) = plot(B{3}(:,XAxisMode),B{3}(:,5),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(B{3}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(B{3}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)),B{3}(z(1),5)-(YAxis(2)-YAxis(1))/30,'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,B{3}(z(1),5)-(YAxis(2)-YAxis(1))/30,'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,B{3}(z(1),5)-(YAxis(2)-YAxis(1))/30,'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 4:size(B,2)
                    B{i}(:,5) = real(asind(Couplant.Velocity/1e3./B{i}(:,4)));
                    b(i) = plot(B{i}(:,XAxisMode),B{i}(:,5),'LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
    else
        if  LambModes
            BLamb{1}(:,5) = real(asind(Couplant.Velocity/1e3./BLamb{1}(:,4)));
            bLamb = plot(BLamb{1}(:,XAxisMode),BLamb{1}(:,5),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(BLamb{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(BLamb{1}(:,1)-(ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1)),BLamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,BLamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel3X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,BLamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            BLamb{2}(:,5) = real(asind(Couplant.Velocity/1e3./BLamb{2}(:,4)));
            bLamb(2) = plot(BLamb{2}(:,XAxisMode),BLamb{2}(:,5),'LineWidth',LineWidth,'Color',BColor);
            if  ModeLabels
                z = find(abs(BLamb{2}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(BLamb{2}(:,1)-(ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                if  XAxisMode == 1
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1)),BLamb{2}(z(1),5)-(YAxis(2)+YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 2
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,BLamb{2}(z(1),5)-(YAxis(2)+YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                elseif XAxisMode == 3
                    text((ModeLabel1X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,BLamb{2}(z(1),5)-(YAxis(2)+YAxis(1))/30,'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                end
            end
            if  HigherOrderModes
                for i = 3:size(BLamb,2)
                    BLamb{i}(:,5) = real(asind(Couplant.Velocity/1e3./BLamb{i}(:,4)));
                    bLamb(i) = plot(BLamb{i}(:,XAxisMode),BLamb{i}(:,5),'LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
        if  ShearHorizontalModes
            if  SuperLayerSize > 1 && ~SymmetricSystem
                BShear{1}(:,5) = real(asind(Couplant.Velocity/1e3./BShear{1}(:,4)));
                bShear = plot(BShear{1}(:,XAxisMode),BShear{1}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                if  ModeLabels
                    z = find(abs(BShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(BShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                    if  XAxisMode == 1
                        text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)),BShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                    elseif XAxisMode == 2
                        text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,BShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                    elseif XAxisMode == 3
                        text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,BShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'B$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                    end
                end
                if  HigherOrderModes
                    for i = 2:size(BShear,2)
                        BShear{i}(:,5) = real(asind(Couplant.Velocity/1e3./BShear{i}(:,4)));
                        bShear(i) = plot(BShear{i}(:,XAxisMode),BShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ShearHorizontalModes && ~isempty(SShear{1})
                    SShear{1}(:,5) = real(asind(Couplant.Velocity/1e3./SShear{1}(:,4)));
                    sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                    if  ModeLabels
                        z = find(abs(SShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(SShear{1}(:,1)-(ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)))));
                        if  XAxisMode == 1
                            text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1)),SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                        elseif XAxisMode == 2
                            text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                        elseif XAxisMode == 3
                            text((ModeLabel2X*(XAxis(2)-XAxis(1))+XAxis(1))*PlateThickness,SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
                        end
                    end
                    if  HigherOrderModes
                        for i = 2:size(SShear,2)
                            SShear{i}(:,5) = real(asind(Couplant.Velocity/1e3./SShear{i}(:,4)));
                            sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                        end
                    end
                end
                if  AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                    for i = 1:size(AShear,2)
                        AShear{i}(:,5) = real(asind(Couplant.Velocity/1e3./AShear{i}(:,4)));
                        aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
                    end
                end    
            end
        end
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            BScholte{i}(:,5) = real(asind(Couplant.Velocity/1e3./BScholte{i}(:,4)));
            bScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
        end
    end    
end
ax = gca;
ax.Box = 'on';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 1
    if  XAxisMode ~= 3
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_'))];
    else
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_'))];
    end
elseif ~SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        end
    end
end
if  HeadLine > 0
    if  FluidLoading == 1
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
        end
    end
    ax.Title.String = String;
end
if  XAxisMode == 1
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = XAxis;
elseif XAxisMode == 2
    ax.XLabel.String = 'Frequency (MHz)';
    ax.XLim = XAxis/1e3;
elseif XAxisMode == 3
    ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
    ax.XLim = XAxis*PlateThickness;
end
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = ['Coincidence angle in ',char(join(split(Couplant.Name,'_'),'\_')),' ($^{\circ}$)'];
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLim = YAxis;
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            if  ~Crop
                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[50 30])
                print(f,fullfile(Directory,FileName),'-dpdf')
            else
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
            end
        end
        if  PNG
            if  ~Crop
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
            else
                exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
            end
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
        return
    end
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
% c1 = d.UIContextMenu;
function output_txt = Cursor(~,event_obj)
    if  Symmetric
        if  ~Decoupled
            if  ~isempty(S{1}) && all(event_obj.Target.Color == SColor)
                for i = 1:length(S)
                    if  i == 1 && event_obj.Target.YData(1) == S{1}(1,5)
                        ModeName = 'S$_0$';
                        break
                    elseif i == 2 && event_obj.Target.YData(1) == S{2}(1,5)
                        ModeName = 'S$_1$';
                        break
                    elseif i > 2 && (event_obj.Target.XData(1) == S{i}(1,1) || event_obj.Target.XData(1) == S{i}(1,2) || event_obj.Target.XData(1) == S{i}(1,3))
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(A{1}) && all(event_obj.Target.Color == AColor)
                for i = 1:length(A)
                    if  event_obj.Target.XData(1) == A{i}(1,1) || event_obj.Target.XData(1) == A{i}(1,2) || event_obj.Target.XData(1) == A{i}(1,3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
        else
            if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(SLamb)
                    if  event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3)
                        ModeName = ['S$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(SShear)
                    if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                        ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(ALamb)
                    if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                        ModeName = ['A$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                for i = 1:length(AShear)
                    if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                        ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                        break
                    end
                end
            end
        end
        if  ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(SScholte)
                if  event_obj.Target.XData(1) == SScholte{i}(1,1) || event_obj.Target.XData(1) == SScholte{i}(1,2) || event_obj.Target.XData(1) == SScholte{i}(1,3)
                    ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
        if  ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(AScholte)
                if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                    ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
    else
        if  ~Decoupled
            if  ~isempty(B{1})
                for i = 1:length(B)
                    if  i == 1 && event_obj.Target.YData(1) == B{1}(1,5)
                        ModeName = 'B$_0$';
                        break 
                    elseif i == 2 && event_obj.Target.YData(1) == B{2}(1,5)
                        ModeName = 'B$_1$';
                        break    
                    elseif i == 3 && event_obj.Target.YData(1) == B{3}(1,5)
                        ModeName = 'B$_2$';
                        break
                    elseif i > 3 && (event_obj.Target.XData(1) == B{i}(1,1) || event_obj.Target.XData(1) == B{i}(1,2) || event_obj.Target.XData(1) == B{i}(1,3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end                
        else
            if  ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
                for i = 1:length(BLamb)
                    if  i == 1 && event_obj.Target.YData(1) == BLamb{1}(1,5)
                        ModeName = 'B$_0$';
                        break
                    elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,5)
                        ModeName = 'B$_1$';
                        break
                    elseif i > 2 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                        ModeName = ['B$_{',num2str(i-1),'}$'];
                        break
                    end
                end
            end
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  ~isempty(BShear{1}) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(BShear)
                        if  event_obj.Target.XData(1) == BShear{i}(1,1) || event_obj.Target.XData(1) == BShear{i}(1,2) || event_obj.Target.XData(1) == BShear{i}(1,3)
                            ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            break
                        end
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(SShear)
                        if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                            ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                            break
                        end
                    end
                end
                if  ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
                    for i = 1:length(AShear)
                        if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                            ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                            break
                        end
                    end
                end
            end
        end
        if  ~isempty(BScholte{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(BScholte)
                if  i == 1 && event_obj.Target.YData(1) == BScholte{1}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_0$';
                    break
                elseif i == 2 && event_obj.Target.YData(1) == BScholte{2}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_1$';
                    break
                elseif i > 2 && (event_obj.Target.XData(1) == BScholte{i}(1,1) || event_obj.Target.XData(1) == BScholte{i}(1,2) || event_obj.Target.XData(1) == BScholte{i}(1,3))
                    ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end        
    end
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$\theta_{\mathrm{C}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$\theta_{\mathrm{C}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f\cdot d$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$\theta_{\mathrm{C}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
%     c2 = uicontextmenu;
%     d.UIContextMenu = c2;
%     uimenu(c2,'Label','Hide','Callback',@State_Callback);
%     function State_Callback(~,~)
%         event_obj.Target.LineStyle  = 'none';
%         d.UIContextMenu = c1;
%     end        
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 60],'CloseRequestFcn',@CloseRequest);
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  XAxisMode == 1 % (kHz)
        Value = .1*XAxis(2);
        String = {'Frequency (kHz)',['Coincidence angle (',char(176),')']};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)',['Coincidence angle (',char(176),')']};
    elseif XAxisMode == 3 % (MHz*mm)
        Value = .1*XAxis(2)*PlateThickness;
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],['Coincidence angle (',char(176),')']};
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'TooltipString','Select constant quantity.','String',String,'Position',[10 20 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'TooltipString','Enter a value for the above selected quantity.','Position',[150 20 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                Value = .1*XAxis(2)*PlateThickness; 
            end
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % coincidence angle
            Value = .9*YAxis(2);
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function Value_Callback(source,~)
        Value = str2double(source.String);
        if  Mode == 1 % frequency
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        elseif Mode == 2 % coincidence angle
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Symmetric
            if  ~Decoupled
                if  Mode == 1 % frequency
                    if  AntisymmetricModes && ~isempty(A{1})
                        for i = 1:length(a)
                            if  a(i).XData(1) < Value
                                if  a(i).XData(end) >= Value
                                    z = abs(a(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  SymmetricModes && ~isempty(S{1})
                        for i = 1:length(s)
                            if  s(i).XData(1) < Value
                                if  s(i).XData(end) >= Value
                                    z = abs(s(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  aScholte(i).XData(1) < Value
                                if  aScholte(i).XData(end) >= Value
                                    z = abs(aScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  sScholte(i).XData(1) < Value
                                if  sScholte(i).XData(end) >= Value
                                    z = abs(sScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                elseif Mode == 2 % coincidence angle
                    if  AntisymmetricModes && ~isempty(A{1})
                        for i = 1:length(a)
                            if  min(a(i).YData) < Value && max(a(i).YData) > Value
                                z = abs(a(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                            end
                        end
                    end
                    if  SymmetricModes && ~isempty(S{1})
                        for i = 1:length(s)
                            if  min(s(i).YData) < Value && max(s(i).YData) > Value
                                z = abs(s(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  min(aScholte(i).YData) < Value && max(aScholte(i).YData) > Value
                                z = abs(aScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  min(sScholte(i).YData) < Value && max(sScholte(i).YData) > Value
                                z = abs(sScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                            end
                        end
                    end
                end
            else
                if  Mode == 1 % frequency
                    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                        for i = 1:length(aLamb)
                            if  aLamb(i).XData(1) < Value
                                if  aLamb(i).XData(end) >= Value
                                    z = abs(aLamb(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                        for i = 1:length(sLamb)
                            if  sLamb(i).XData(1) < Value
                                if  sLamb(i).XData(end) >= Value
                                    z = abs(sLamb(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                        for i = 1:length(aShear)
                            if  aShear(i).XData(1) < Value
                                if  aShear(i).XData(end) >= Value
                                    z = abs(aShear(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                        for i = 1:length(sShear)
                            if  sShear(i).XData(1) < Value
                                if  sShear(i).XData(end) >= Value
                                    z = abs(sShear(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  aScholte(i).XData(1) < Value
                                if  aScholte(i).XData(end) >= Value
                                    z = abs(aScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  sScholte(i).XData(1) < Value
                                if  sScholte(i).XData(end) >= Value
                                    z = abs(sScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                elseif Mode == 2 % coincidence angle
                    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                        for i = 1:length(aLamb)
                            if  min(aLamb(i).YData) < Value && max(aLamb(i).YData) > Value
                                z = abs(aLamb(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                            end
                        end
                    end
                    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                        for i = 1:length(sLamb)
                            if  min(sLamb(i).YData) < Value && max(sLamb(i).YData) > Value
                                z = abs(sLamb(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                            end
                        end
                    end
                    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                        for i = 1:length(aShear)
                            if  min(aShear(i).YData) < Value && max(aShear(i).YData) > Value
                                z = abs(aShear(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                            end
                        end
                    end
                    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                        for i = 1:length(sShear)
                            if  min(sShear(i).YData) < Value && max(sShear(i).YData) > Value
                                z = abs(sShear(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                        for i = 1:length(aScholte)
                            if  min(aScholte(i).YData) < Value && max(aScholte(i).YData) > Value
                                z = abs(aScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                        for i = 1:length(sScholte)
                            if  min(sScholte(i).YData) < Value && max(sScholte(i).YData) > Value
                                z = abs(sScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                            end
                        end
                    end
                end
            end
        else
            if  ~Decoupled
                if  Mode == 1 % frequency
                    if  ~isempty(B{1})
                        for i = 1:length(b)
                            if  b(i).XData(1) < Value
                                if  b(i).XData(end) >= Value
                                    z = abs(b(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end 
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  bScholte(i).XData(1) < Value
                                if  bScholte(i).XData(end) >= Value
                                    z = abs(bScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                elseif Mode == 2 % coincidence angle
                    if  ~isempty(B{1})
                        for i = 1:length(b)
                            if  min(b(i).YData) < Value && max(b(i).YData) > Value
                                z = abs(b(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  min(bScholte(i).YData) < Value && max(bScholte(i).YData) > Value
                                z = abs(bScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                            end
                        end
                    end
                end
            else
                if  Mode == 1 % frequency
                    if  LambModes && ~isempty(BLamb{1})
                        for i = 1:length(bLamb)
                            if  bLamb(i).XData(1) < Value
                                if  bLamb(i).XData(end) >= Value
                                    z = abs(bLamb(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bLamb(i),bLamb(i).XData(q),bLamb(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                    if  SuperLayerSize > 1 && ~SymmetricSystem
                        if  ShearHorizontalModes && ~isempty(BShear{1})
                            for i = 1:length(bShear)
                                if  bShear(i).XData(1) < Value
                                    if  bShear(i).XData(end) >= Value
                                        z = abs(bShear(i).XData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(bShear(i),bShear(i).XData(q),bShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                    elseif SuperLayerSize == 1 || SymmetricSystem
                        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                            for i = 1:length(aShear)
                                if  aShear(i).XData(1) < Value
                                    if  aShear(i).XData(end) >= Value
                                        z = abs(aShear(i).XData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                            for i = 1:length(sShear)
                                if  sShear(i).XData(1) < Value
                                    if  sShear(i).XData(end) >= Value
                                        z = abs(sShear(i).XData-Value);
                                        q = find(z == min(z),1);
                                        dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                    end
                                else
                                    break
                                end
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  bScholte(i).XData(1) < Value
                                if  bScholte(i).XData(end) >= Value
                                    z = abs(bScholte(i).XData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                                end
                            else
                                break
                            end
                        end
                    end
                elseif Mode == 2 % coincidence angle
                    if  LambModes && ~isempty(BLamb{1})
                        for i = 1:length(bLamb)
                            if  min(bLamb(i).YData) < Value && max(bLamb(i).YData) > Value
                                z = abs(bLamb(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bLamb(i),bLamb(i).XData(q),bLamb(i).YData(q));
                            end
                        end
                    end
                    if  SuperLayerSize > 1 && ~SymmetricSystem
                        if  ShearHorizontalModes && ~isempty(BShear{1})
                            for i = 1:length(bShear)
                                if  min(bShear(i).YData) < Value && max(bShear(i).YData) > Value
                                    z = abs(bShear(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(bShear(i),bShear(i).XData(q),bShear(i).YData(q));
                                end
                            end
                        end
                    elseif SuperLayerSize == 1 || SymmetricSystem
                        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                            for i = 1:length(aShear)
                                if  min(aShear(i).YData) < Value && max(aShear(i).YData) > Value
                                    z = abs(aShear(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                                end
                            end
                        end
                        if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                            for i = 1:length(sShear)
                                if  min(sShear(i).YData) < Value && max(sShear(i).YData) > Value
                                    z = abs(sShear(i).YData-Value);
                                    q = find(z == min(z),1);
                                    dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                                end
                            end
                        end
                    end
                    if  ScholteModes && ~isempty(BScholte{1})
                        for i = 1:length(bScholte)
                            if  min(bScholte(i).YData) < Value && max(bScholte(i).YData) > Value
                                z = abs(bScholte(i).YData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                            end
                        end
                    end
                end
            end
        end
    end
    function CloseRequest(~,~)
        delete(l)
        delete(dt)
        delete(f_Analyze)
    end
end
end