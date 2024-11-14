% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% =========================================================================
function Export_Polar(a,XLSX,TXT,MAT)
%#ok<*AGROW>
FrequencyRange = 0:a.FrequencyResolution_Polar:a.FrequencyLimit_Polar;
FrequencyRange(1) = a.FrequencyRangeStart2;
q = find(FrequencyRange == a.Frequency_Polar);
if  isempty(q)
    if  a.Frequency_Polar > ceil(max(FrequencyRange)) || a.Frequency_Polar < 0
        errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(ceil(max(FrequencyRange))),' kHz.'],'Error');
        return
    else
        [~,q] = min(abs(FrequencyRange-a.Frequency_Polar));
        a.Frequency_Polar = FrequencyRange(q);
    end
end
if  a.PropagationAngleMode_Polar == 1
    Limit = 180;
elseif a.PropagationAngleMode_Polar == 2
    Limit = 90;
end
PropagationAngle = 0:a.PropagationAngleStep_Polar:Limit;
for n = 1:Limit/a.PropagationAngleStep_Polar+1
    if  a.A0_Polar
        A0(n,1) = a.A_Polar{n,1}(q,1)/1e3; % cp
        A0(n,2) = a.A_Polar{n,1}(q,2)/1e3; % cg1
        A0(n,4) = a.Distance_Polar./a.A_Polar{n,1}(q,2)*1e3; % t1
        A0(n,5) = real(asind(a.Couplant_Polar.Velocity./a.A_Polar{n,1}(q,1))); % excitation angle
        A0(n,6) = a.A_Polar{n,1}(q,1)./a.Frequency_Polar; % wavelength
        A0(n,7) = 2*pi*a.Frequency_Polar/a.A_Polar{n,1}(q,1); % wavenumber
        SkewAngle(n,1) = -atand(a.A_Polar{n,1}(q,3)/a.A_Polar{n,1}(q,2)); % skew angle
        A0_2(n,2) = sqrt(a.A_Polar{n,1}(q,2)^2+a.A_Polar{n,1}(q,3)^2)/1e3; % cAbs
        A0_2(n,3) = a.Distance_Polar/sqrt(a.A_Polar{n,1}(q,2)^2+a.A_Polar{n,1}(q,3)^2)*1e3; % tAbs
    end
    if  a.SH0_Polar
        SH0(n,1) = a.A_Polar{n,2}(q,1)/1e3;
        SH0(n,2) = a.A_Polar{n,2}(q,2)/1e3;
        SH0(n,4) = a.Distance_Polar./a.A_Polar{n,2}(q,2)*1e3;
        SH0(n,5) = real(asind(a.Couplant_Polar.Velocity./a.A_Polar{n,2}(q,1)));
        SH0(n,6) = a.A_Polar{n,2}(q,1)./a.Frequency_Polar;
        SH0(n,7) = 2*pi*a.Frequency_Polar/a.A_Polar{n,2}(q,1);
        SkewAngle(n,2) = -atand(a.A_Polar{n,2}(q,3)/a.A_Polar{n,2}(q,2));
        SH0_2(n,2) = sqrt(a.A_Polar{n,2}(q,2)^2+a.A_Polar{n,2}(q,3)^2)/1e3;
        SH0_2(n,3) = a.Distance_Polar/sqrt(a.A_Polar{n,2}(q,2)^2+a.A_Polar{n,2}(q,3)^2)*1e3;
    end
    if  a.S0_Polar
        S0(n,1) = a.A_Polar{n,3}(q,1)/1e3;
        S0(n,2) = a.A_Polar{n,3}(q,2)/1e3;
        S0(n,4) = a.Distance_Polar./a.A_Polar{n,3}(q,2)*1e3;
        S0(n,5) = real(asind(a.Couplant_Polar.Velocity./a.A_Polar{n,3}(q,1)));
        S0(n,6) = a.A_Polar{n,3}(q,1)./a.Frequency_Polar;
        S0(n,7) = 2*pi*a.Frequency_Polar/a.A_Polar{n,3}(q,1);
        SkewAngle(n,3) = -atand(a.A_Polar{n,3}(q,3)/a.A_Polar{n,3}(q,2));
        S0_2(n,2) = sqrt(a.A_Polar{n,3}(q,2)^2+a.A_Polar{n,3}(q,3)^2)/1e3;
        S0_2(n,3) = a.Distance_Polar/sqrt(a.A_Polar{n,3}(q,2)^2+a.A_Polar{n,3}(q,3)^2)*1e3;
    end
end
if  abs(SkewAngle(1,1)) < .5
    SkewAngle(1,:) = 0;
    SkewAngle(end,:) = 0;
    if  a.PropagationAngleMode_Polar == 1
        SkewAngle(ceil(.5*length(SkewAngle)),:) = 0;
    end
end
if  a.A0_Polar 
    if  a.PropagationAngleMode_Polar == 1
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,1));
        A0 = vertcat(A0,A0(2:end,:));
        A0(:,3) = vertcat(SkewAngle(:,1),SkewAngle(2:end,1));
        A0_2 = vertcat(A0_2,A0_2(2:end,:));
        A0_2(:,1) = vertcat(RayAngle,pi+RayAngle(2:end));
    elseif a.PropagationAngleMode_Polar == 2
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,1));
        A0 = vertcat(A0,flipud(A0(2:end-1,:)),A0(1:end,:),flipud(A0(1:end-1,:)));
        A0(:,3) = vertcat(SkewAngle(:,1),-flipud(SkewAngle(2:end-1,1)),SkewAngle(1:end,1),-flipud(SkewAngle(1:end-1,1)));
        A0_2 = vertcat(A0_2,flipud(A0_2(2:end-1,:)),A0_2(1:end,:),flipud(A0_2(1:end-1,:)));
        A0_2(:,1) = vertcat(RayAngle,pi-flipud(RayAngle(1:end-1)),pi+RayAngle(2:end),2*pi-flipud(RayAngle(1:end-1)));
    end
end
if  a.SH0_Polar
    if  a.PropagationAngleMode_Polar == 1
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,2));
        SH0 = vertcat(SH0,SH0(2:end,:));
        SH0(:,3) = vertcat(SkewAngle(:,2),SkewAngle(2:end,2));
        SH0_2 = vertcat(SH0_2,SH0_2(2:end,:));
        SH0_2(:,1) = vertcat(RayAngle,pi+RayAngle(2:end));
    elseif a.PropagationAngleMode_Polar == 2
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,2));
        SH0 = vertcat(SH0,flipud(SH0(2:end-1,:)),SH0(1:end,:),flipud(SH0(1:end-1,:)));
        SH0(:,3) = vertcat(SkewAngle(:,2),-flipud(SkewAngle(2:end-1,2)),SkewAngle(1:end,2),-flipud(SkewAngle(1:end-1,2)));
        SH0_2 = vertcat(SH0_2,flipud(SH0_2(2:end-1,:)),SH0_2(1:end,:),flipud(SH0_2(1:end-1,:)));
        SH0_2(:,1) = vertcat(RayAngle,pi-flipud(RayAngle(1:end-1)),pi+RayAngle(2:end),2*pi-flipud(RayAngle(1:end-1)));
    end
end
if  a.S0_Polar
    if  a.PropagationAngleMode_Polar == 1
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,3));
        S0 = vertcat(S0,S0(2:end,:));
        S0(:,3) = vertcat(SkewAngle(:,3),SkewAngle(2:end,3));
        S0_2 = vertcat(S0_2,S0_2(2:end,:));
        S0_2(:,1) = vertcat(RayAngle,pi+RayAngle(2:end));
    elseif a.PropagationAngleMode_Polar == 2
        RayAngle = deg2rad(PropagationAngle'-SkewAngle(:,3));
        S0 = vertcat(S0,flipud(S0(2:end-1,:)),S0(1:end,:),flipud(S0(1:end-1,:)));
        S0(:,3) = vertcat(SkewAngle(:,3),-flipud(SkewAngle(2:end-1,3)),SkewAngle(1:end,3),-flipud(SkewAngle(1:end-1,3)));
        S0_2 = vertcat(S0_2,flipud(S0_2(2:end-1,:)),S0_2(1:end,:),flipud(S0_2(1:end-1,:)));
        S0_2(:,1) = vertcat(RayAngle,pi-flipud(RayAngle(1:end-1)),pi+RayAngle(2:end),2*pi-flipud(RayAngle(1:end-1)));
    end
end
if  a.PropagationAngleMode_Polar == 1
    z = 0:2*pi/(2*(n-1)):2*pi;
elseif a.PropagationAngleMode_Polar == 2
    z = 0:2*pi/(4*(n-1)):2*pi;
end
if  a.A0_Polar && a.SH0_Polar && a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','A0 Phase velocity (m/ms)','A0 Energy velocity1 (m/ms)','A0 Skew angle (deg)','A0 Propagation time1 (micsec)','A0 Coincidence angle (deg)','A0 Wavelength (mm)','A0 Wavenumber (rad/mm)','S0 Phase velocity (m/ms)','S0 Energy velocity1 (m/ms)','S0 Skew angle (deg)','S0 Propagation time1 (micsec)','S0 Coincidence angle (deg)','S0 Wavelength (mm)','S0 Wavenumber (rad/mm)','S1 Phase velocity (m/ms)','S1 Energy velocity1 (m/ms)','S1 Skew angle (deg)','S1 Propagation time1 (micsec)','S1 Coincidence angle (deg)','S1 Wavelength (mm)','S1 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'A0 Ray angle (rad)','A0 Energy velocity absolute (m/ms)','A0 Propagation time absolute (micsec)','S0 Ray angle (rad)','S0 Energy velocity absolute (m/ms)','S0 Propagation time absolute (micsec)','S1 Ray angle (rad)','S1 Energy velocity absolute (m/ms)','S1 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','B0 Phase velocity (m/ms)','B0 Energy velocity1 (m/ms)','B0 Skew angle (deg)','B0 Propagation time1 (micsec)','B0 Coincidence angle (deg)','B0 Wavelength (mm)','B0 Wavenumber (rad/mm)','B1 Phase velocity (m/ms)','B1 Energy velocity1 (m/ms)','B1 Skew angle (deg)','B1 Propagation time1 (micsec)','B1 Coincidence angle (deg)','B1 Wavelength (mm)','B1 Wavenumber (rad/mm)','B2 Phase velocity (m/ms)','B2 Energy velocity1 (m/ms)','B2 Skew angle (deg)','B2 Propagation time1 (micsec)','B2 Coincidence angle (deg)','B2 Wavelength (mm)','B2 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'B0 Ray angle (rad)','B0 Energy velocity absolute (m/ms)','B0 Propagation time absolute (micsec)','B1 Ray angle (rad)','B1 Energy velocity absolute (m/ms)','B1 Propagation time absolute (micsec)','B2 Ray angle (rad)','B2 Energy velocity absolute (m/ms)','B2 Propagation time absolute (micsec)'});
    end
elseif a.A0_Polar && a.SH0_Polar && ~a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),'VariableNames',{'Propagation angle (rad)','A0 Phase velocity (m/ms)','A0 Energy velocity1 (m/ms)','A0 Skew angle (deg)','A0 Propagation time1 (micsec)','A0 Coincidence angle (deg)','A0 Wavelength (mm)','A0 Wavenumber (rad/mm)','S0 Phase velocity (m/ms)','S0 Energy velocity1 (m/ms)','S0 Skew angle (deg)','S0 Propagation time1 (micsec)','S0 Coincidence angle (deg)','S0 Wavelength (mm)','S0 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),'VariableNames',{'A0 Ray angle (rad)','A0 Energy velocity absolute (m/ms)','A0 Propagation time absolute (micsec)','S0 Ray angle (rad)','S0 Energy velocity absolute (m/ms)','S0 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),'VariableNames',{'Propagation angle (rad)','B0 Phase velocity (m/ms)','B0 Energy velocity1 (m/ms)','B0 Skew angle (deg)','B0 Propagation time1 (micsec)','B0 Coincidence angle (deg)','B0 Wavelength (mm)','B0 Wavenumber (rad/mm)','B1 Phase velocity (m/ms)','B1 Energy velocity1 (m/ms)','B1 Skew angle (deg)','B1 Propagation time1 (micsec)','B1 Coincidence angle (deg)','B1 Wavelength (mm)','B1 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),'VariableNames',{'B0 Ray angle (rad)','B0 Energy velocity absolute (m/ms)','B0 Propagation time absolute (micsec)','B1 Ray angle (rad)','B1 Energy velocity absolute (m/ms)','B1 Propagation time absolute (micsec)'});
    end
elseif a.A0_Polar && ~a.SH0_Polar && a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','A0 Phase velocity (m/ms)','A0 Energy velocity1 (m/ms)','A0 Skew angle (deg)','A0 Propagation time1 (micsec)','A0 Coincidence angle (deg)','A0 Wavelength (mm)','A0 Wavenumber (rad/mm)','S1 Phase velocity (m/ms)','S1 Energy velocity1 (m/ms)','S1 Skew angle (deg)','S1 Propagation time1 (micsec)','S1 Coincidence angle (deg)','S1 Wavelength (mm)','S1 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'A0 Ray angle (rad)','A0 Energy velocity absolute (m/ms)','A0 Propagation time absolute (micsec)','S1 Ray angle (rad)','S1 Energy velocity absolute (m/ms)','S1 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','B0 Phase velocity (m/ms)','B0 Energy velocity1 (m/ms)','B0 Skew angle (deg)','B0 Propagation time1 (micsec)','B0 Coincidence angle (deg)','B0 Wavelength (mm)','B0 Wavenumber (rad/mm)','B2 Phase velocity (m/ms)','B2 Energy velocity1 (m/ms)','B2 Skew angle (deg)','B2 Propagation time1 (micsec)','B2 Coincidence angle (deg)','B2 Wavelength (mm)','B2 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'B0 Ray angle (rad)','B0 Energy velocity absolute (m/ms)','B0 Propagation time absolute (micsec)','B2 Ray angle (rad)','B2 Energy velocity absolute (m/ms)','B2 Propagation time absolute (micsec)'});
    end
elseif ~a.A0_Polar && a.SH0_Polar && a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','S0 Phase velocity (m/ms)','S0 Energy velocity1 (m/ms)','S0 Skew angle (deg)','S0 Propagation time1 (micsec)','S0 Coincidence angle (deg)','S0 Wavelength (mm)','S0 Wavenumber (rad/mm)','S1 Phase velocity (m/ms)','S1 Energy velocity1 (m/ms)','S1 Skew angle (deg)','S1 Propagation time1 (micsec)','S1 Coincidence angle (deg)','S1 Wavelength (mm)','S1 Wavenumber (rad/mm)'});
        Table2 = table(SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'S0 Ray angle (rad)','S0 Energy velocity absolute (m/ms)','S0 Propagation time absolute (micsec)','S1 Ray angle (rad)','S1 Energy velocity absolute (m/ms)','S1 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','B1 Phase velocity (m/ms)','B1 Energy velocity1 (m/ms)','B1 Skew angle (deg)','B1 Propagation time1 (micsec)','B1 Coincidence angle (deg)','B1 Wavelength (mm)','B1 Wavenumber (rad/mm)','B2 Phase velocity (m/ms)','B2 Energy velocity1 (m/ms)','B2 Skew angle (deg)','B2 Propagation time1 (micsec)','B2 Coincidence angle (deg)','B2 Wavelength (mm)','B2 Wavenumber (rad/mm)'});
        Table2 = table(SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'B1 Ray angle (rad)','B1 Energy velocity absolute (m/ms)','B1 Propagation time absolute (micsec)','B2 Ray angle (rad)','B2 Energy velocity absolute (m/ms)','B2 Propagation time absolute (micsec)'});
    end
elseif a.A0_Polar && ~a.SH0_Polar && ~a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),'VariableNames',{'Propagation angle (rad)','A0 Phase velocity (m/ms)','A0 Energy velocity1 (m/ms)','A0 Skew angle (deg)','A0 Propagation time1 (micsec)','A0 Coincidence angle (deg)','A0 Wavelength (mm)','A0 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),'VariableNames',{'A0 Ray angle (rad)','A0 Energy velocity absolute (m/ms)','A0 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',A0(:,1),A0(:,2),A0(:,3),A0(:,4),A0(:,5),A0(:,6),A0(:,7),'VariableNames',{'Propagation angle (rad)','B0 Phase velocity (m/ms)','B0 Energy velocity1 (m/ms)','B0 Skew angle (deg)','B0 Propagation time1 (micsec)','B0 Coincidence angle (deg)','B0 Wavelength (mm)','B0 Wavenumber (rad/mm)'});
        Table2 = table(A0_2(:,1),A0_2(:,2),A0_2(:,3),'VariableNames',{'B0 Ray angle (rad)','B0 Energy velocity absolute (m/ms)','B0 Propagation time absolute (micsec)'});
    end
elseif ~a.A0_Polar && a.SH0_Polar && ~a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),'VariableNames',{'Propagation angle (rad)','S0 Phase velocity (m/ms)','S0 Energy velocity1 (m/ms)','S0 Skew angle (deg)','S0 Propagation time1 (micsec)','S0 Coincidence angle (deg)','S0 Wavelength (mm)','S0 Wavenumber (rad/mm)'});
        Table2 = table(SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),'VariableNames',{'S0 Ray angle (rad)','S0 Energy velocity absolute (m/ms)','S0 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',SH0(:,1),SH0(:,2),SH0(:,3),SH0(:,4),SH0(:,5),SH0(:,6),SH0(:,7),'VariableNames',{'Propagation angle (rad)','B1 Phase velocity (m/ms)','B1 Energy velocity1 (m/ms)','B1 Skew angle (deg)','B1 Propagation time1 (micsec)','B1 Coincidence angle (deg)','B1 Wavelength (mm)','B1 Wavenumber (rad/mm)'});
        Table2 = table(SH0_2(:,1),SH0_2(:,2),SH0_2(:,3),'VariableNames',{'B1 Ray angle (rad)','B1 Energy velocity absolute (m/ms)','B1 Propagation time absolute (micsec)'});
    end
elseif ~a.A0_Polar && ~a.SH0_Polar && a.S0_Polar
    if  isscalar(a.LayerOrientations_Polar) || a.SymmetricSystem_Polar
        Table1 = table(z',S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','S1 Phase velocity (m/ms)','S1 Energy velocity1 (m/ms)','S1 Skew angle (deg)','S1 Propagation time1 (micsec)','S1 Coincidence angle (deg)','S1 Wavelength (mm)','S1 Wavenumber (rad/mm)'});
        Table2 = table(S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'S1 Ray angle (rad)','S1 Energy velocity absolute (m/ms)','S1 Propagation time absolute (micsec)'});
    else
        Table1 = table(z',S0(:,1),S0(:,2),S0(:,3),S0(:,4),S0(:,5),S0(:,6),S0(:,7),'VariableNames',{'Propagation angle (rad)','B2 Phase velocity (m/ms)','B2 Energy velocity1 (m/ms)','B2 Skew angle (deg)','B2 Propagation time1 (micsec)','B2 Coincidence angle (deg)','B2 Wavelength (mm)','B2 Wavenumber (rad/mm)'});
        Table2 = table(S0_2(:,1),S0_2(:,2),S0_2(:,3),'VariableNames',{'B2 Ray angle (rad)','B2 Energy velocity absolute (m/ms)','B2 Propagation time absolute (micsec)'});
    end
end
try
    if  XLSX
        writetable(Table1,fullfile(a.Directory4,[a.FileName_Polar,'_Profiles@',num2str(a.Frequency_Polar*a.PlateThickness_Polar/1e3),'MHzmm_A.xlsx']))
        writetable(Table2,fullfile(a.Directory4,[a.FileName_Polar,'_Profiles@',num2str(a.Frequency_Polar*a.PlateThickness_Polar/1e3),'MHzmm_B.xlsx']))
    end
    if  TXT
        writetable(Table1,fullfile(a.Directory4,[a.FileName_Polar,'_Profiles@',num2str(a.Frequency_Polar*a.PlateThickness_Polar/1e3),'MHzmm_A.txt']))
        writetable(Table2,fullfile(a.Directory4,[a.FileName_Polar,'_Profiles@',num2str(a.Frequency_Polar*a.PlateThickness_Polar/1e3),'MHzmm_B.txt']))
    end
    if  MAT
        fd = a.Frequency_Polar*a.PlateThickness_Polar/1e3;
        if  fd > 1e-4
            M = matfile(fullfile(a.Directory4,[a.FileName_Polar,'_Profiles_',replace(sprintf('%g',fd),'.','p'),'MHzmm']),'Writable',true); %#ok<*NASGU>
            eval(sprintf('M.Profiles_%sMHzmm_A = Table1;',replace(sprintf('%g',fd),'.','p')))
            eval(sprintf('M.Profiles_%sMHzmm_B = Table2;',replace(sprintf('%g',fd),'.','p')))
        else
            M = matfile(fullfile(a.Directory4,[a.FileName_Polar,'_Profiles_',replace(sprintf('%f',fd),'.','p'),'MHzmm']),'Writable',true);
            eval(sprintf('M.Profiles_%sMHzmm_A = Table1;',replace(sprintf('%f',fd),'.','p')))
            eval(sprintf('M.Profiles_%sMHzmm_B = Table2;',replace(sprintf('%f',fd),'.','p')))
        end
    end
catch ME
    st = dbstack;
    level = find(matches({ME.stack.name},st(1).name));
    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
    return
end