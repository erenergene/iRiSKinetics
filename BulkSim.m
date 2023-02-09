clc;clear;close all

addpath(genpath(pwd))
load("SimOutputs/RefSimData.mat") % Loads variables needed for simulations

%% Get R_Si and R_T1 for air

R_Si_air_red = Refmat_air(find(lambda==cw_r),find(L==0)); % Reflectance at Si at air (n = 1.33)
R_T1_air_red = Refmat_air(find(lambda==cw_r),find(L==100):find(L==120)); % Reflectance at T1 = 100 nm to 120 nm oxide at air (n = 1.33)
R_Si_air_blue = Refmat_air(find(lambda==cw_b),find(L==0)); % Reflectance at Si at air (n = 1.33)
R_T1_air_blue = Refmat_air(find(lambda==cw_b),find(L==100):find(L==120)); % Reflectance at T1 = 100 nm to 120 nm oxide at air (n = 1.33)
R_Si_air_green = Refmat_air(find(lambda==cw_g),find(L==0)); % Reflectance at Si at air (n = 1.33)
R_T1_air_green = Refmat_air(find(lambda==cw_g),find(L==100):find(L==120)); % Reflectance at T1 = 100 nm to 120 nm oxide at air (n = 1.33)

%% Get ratio R_Si - R_T1 / R_Si + R_T1

subtrair_red = R_Si_air_red - R_T1_air_red; % Subtract T1 oxide reflectance from Silicon reflectance
addair_red = R_Si_air_red + R_T1_air_red; % Add T1 oxide reflectance and Silicon reflectance
rat_red = subtrair_red ./ addair_red; % Ratio 
subtrair_blue = R_Si_air_blue - R_T1_air_blue; % Subtract T1 oxide reflectance from Silicon reflectance
addair_blue = R_Si_air_blue + R_T1_air_blue; % Add T1 oxide reflectance and Silicon reflectance
rat_blue = subtrair_blue ./ addair_blue; % Ratio 
subtrair_green = R_Si_air_green - R_T1_air_green; % Subtract T1 oxide reflectance from Silicon reflectance
addair_green = R_Si_air_green + R_T1_air_green; % Add T1 oxide reflectance and Silicon reflectance
rat_green = subtrair_green ./ addair_green; % Ratio 
%% Display Ratio as a function of thickness

figure(1)
hold on
plot(linspace(100,120,201),rat_red,'r','LineWidth',2)
plot(linspace(100,120,201),rat_green,'g','LineWidth',2)
plot(linspace(100,120,201),rat_blue,'b','LineWidth',2)
xlabel('Thickness T1 (nm)'); ylabel('Ratio')
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (FOR AIR)')
saveas(figure(1),[pwd '/Figures/BulkSim/1RatAir.fig']);


%% Define n and L vars

Bulkn = linspace(1.33,1.35,21); % Refractive index from n = 1.33 to 1.35
BulkL = linspace(100,120,201); % Thickness from 100 nm to 120 nm

%% Get reflectance matrices

%% Get reflectance vals: lambda = 400 - 700 nm, thickness = 100 - 120 nm, ref index = 1.33 - 1.35 

%%
% Refmat_lam_L_n_1 = []
% Z1_lam_L_n_1 = []
% 
% tic
% 
% for i = 1:numel(lambda)
%     for j = 1:numel(BulkL)
%         for k = 1:numel(Bulkn)
%             fprintf("Now running lambda %.0f thickness %.0f ref index %.0f\n",i,j,k)
%             [Refmat_lam_L_n_1(i,j,k),Z1_lam_L_n_1(i,j,k)] = multidiel1([Bulkn(k);n_SiO2(i,2);n_Si(i,2)],BulkL(j).*n_SiO2(i,2),lambda(i));
%         end
%     end
% end
% 
% toc
%
% Refmat_lam_L_n = conj(Refmat_lam_L_n_1).*Refmat_lam_L_n_1;

load("Refmat_lam_L_n.mat") % For Faster Code

%% Get reflectance vals: lambda = 460 nm (cw_b) & 633 nm (cw_r), thickness = 0 nm, ref index = 1.33 - 1.35 

Refmat_lamR_L0_n = []
Refmat_lamB_L0_n = []

for i = 1
    for j = 1
        for k = 1:numel(Bulkn)
            [Refmat1_lamR_L0_n(i,j,k),Z1_lamR_L0_n(i,j,k)] = multidiel1([Bulkn(k);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
            [Refmat1_lamB_L0_n(i,j,k),Z1_lamB_L0_n(i,j,k)] = multidiel1([Bulkn(k);n_SiO2(find(lambda==cw_b),2);n_Si(find(lambda==cw_b),2)],0,cw_b); 
        end
    end
end

Refmat_lamR_L0_n = conj(Refmat1_lamR_L0_n).*Refmat1_lamR_L0_n;
Refmat_lamB_L0_n = conj(Refmat1_lamB_L0_n).*Refmat1_lamB_L0_n;

%% Display reflectance for 5 different thickness from 100 to 120 nm for red single wavelength lambda = 633 nm


figure(2)
hold on
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_r),find(BulkL==100),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_r),find(BulkL==105),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_r),find(BulkL==110),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_r),find(BulkL==115),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_r),find(BulkL==120),:),[1 numel(Bulkn)]),'LineWidth',2)

legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Reflectance')
xlim([1.33 1.35])
title('Reflectance')
saveas(figure(2),[pwd '/Figures/BulkSim/2RefRvs5L.fig']);

%% Display reflectance for 5 different thickness from 100 to 120 nm for blue single wavelength lambda = 460 nm

figure(3)
hold on
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_b),find(BulkL==100),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_b),find(BulkL==105),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_b),find(BulkL==110),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_b),find(BulkL==115),:),[1 numel(Bulkn)]),'LineWidth',2)
plot(Bulkn,reshape(Refmat_lam_L_n(find(lambda==cw_b),find(BulkL==120),:),[1 numel(Bulkn)]),'LineWidth',2)

legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Reflectance')
xlim([1.33 1.35])
title('Reflectance')
saveas(figure(3),[pwd '/Figures/BulkSim/3RefBvs5L.fig']);

%% Get T0 and T1 - Reshape reflectance matrices

Refmat_lamR_L0_n = reshape(Refmat_lamR_L0_n,[21 1]); % Reflectance at silicon and red sw for n from 1.33 to 1.35
Refmat_lamB_L0_n = reshape(Refmat_lamB_L0_n,[21 1]); % Reflectance at silicon and blue sw for n from 1.33 to 1.35
Refmat_lamR_L_n = reshape(Refmat_lam_L_n(find(lambda == cw_r),:,:),[201 21])'; % Reflectance at oxide and red sw for n 1.33 to 1.35
Refmat_lamB_L_n = reshape(Refmat_lam_L_n(find(lambda == cw_b),:,:),[201 21])'; % Reflectance at oxide and blue sw for n 1.33 to 1.35

%% Get ratio for red sw

subtrSiT1R = Refmat_lamR_L0_n - Refmat_lamR_L_n; %T0 - T1
addSiT1R = Refmat_lamR_L0_n + Refmat_lamR_L_n; % T0 + T1
ratSiT1R = subtrSiT1R ./ addSiT1R; % (T0 - T1) / (T0 + T1)

%% Get ratio for blue sw

subtrSiT1B = Refmat_lamB_L0_n - Refmat_lamB_L_n; %T0 - T1
addSiT1B = Refmat_lamB_L0_n + Refmat_lamB_L_n; % T0 + T1
ratSiT1B = subtrSiT1B ./ addSiT1B; % (T0 - T1) / (T0 + T1)


%% Display ratio (R_Si - R_T1) / (R_Si + R_T1) as a function of n from 1.33 to 1.35

figure(4) 

hold on
plot(Bulkn,ratSiT1R(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratSiT1R(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratSiT1R(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratSiT1R(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratSiT1R(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')
saveas(figure(4),[pwd '/Figures/BulkSim/4RatSiT1R5Lvsn.fig']);
%% 

figure(5) 

hold on
plot(Bulkn,ratSiT1B(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratSiT1B(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratSiT1B(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratSiT1B(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratSiT1B(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (B)')
saveas(figure(5),[pwd '/Figures/BulkSim/5RatSiT1B5Lvsn.fig']);

%% Get and display slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slpratR(i) = (ratSiT1R(end,i) - ratSiT1R(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

for i = 1:numel(BulkL)
    slpratB(i) = (ratSiT1B(end,i) - ratSiT1B(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

%% Display slope of ratio as a function of thickness

figure(6)
hold on
plot(BulkL,slpratR,'LineWidth',2)
plot(BulkL,slpratB,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')
saveas(figure(6),[pwd '/Figures/BulkSim/6SlpRatRBvsL.fig']);

%% Code to visualize n in a broader range (n = 1 to 3)

%% Define broad n

nbroad = linspace(1,3,201); 

%% Get reflectance vals: lambda = 400 - 700 nm, thickness = 100 - 120 nm, ref index = 1 - 3


for i = 1:numel(nbroad)
[Refmat1_lamR_L0_nbroad(i),Z1(i)] = multidiel1([nbroad(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
end

Refmat_lamR_L0_nbroad = conj(Refmat1_lamR_L0_nbroad).*Refmat1_lamR_L0_nbroad;

%% Get reflectance vals from n = 1.33 to 1.35 and thickness 100 120 nm for red led

tic

for i = 1:numel(nbroad)
for j = 1:numel(BulkL)
fprintf("Now running ref index %.0f\n",i)
fprintf("Now running thickness %.0f\n",j)

[Refmat1_lamR_L_nbroad(i,j),Z1(i,j)] = multidiel1([nbroad(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],BulkL(j).*n_SiO2(find(lambda==cw_r),2),cw_r);        

% Calculates reflectance as a function of thickness L from 100 to 120 nm and refractive index n from 1.33 to 1.35 

end
end

toc

Refmat_lamR_L_nbroad = conj(Refmat1_lamR_L_nbroad).*Refmat1_lamR_L_nbroad; %Refmat_Bulk(n,L) 

Refmat_lamR_L_nbroad = Refmat_lamR_L_nbroad';
subbulknbroad = Refmat_lamR_L0_nbroad - Refmat_lamR_L_nbroad; % Subtract T1 oxide reflectance from Silicon reflectance
addbulknbroad = Refmat_lamR_L0_nbroad + Refmat_lamR_L_nbroad; % Add T1 oxide reflectance and Silicon reflectance
ratbulknbroad = (subbulknbroad ./ addbulknbroad)'; % Ratio

%%

figure(6)
hold on
plot(nbroad,ratbulknbroad(:,find(BulkL==100)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==105)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==110)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==115)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(5),[pwd '/Figures/BulkSim/5RatBulktest.fig']);

%% 
%% LED Integrated Response

load("Osram_Spec_Data_edited.mat") % get LEDs
bluespectrum = Spec_DentalBlue_460nm; clear("Spec_DentalBlue_460nm")
greenspectrum = Spec_Green_517nm; clear("Spec_Green_517nm")
redspectrum = Spec_Red_633nm; clear("Spec_Red_633nm")

intbluespectrum = interp1(Osram_lambda,bluespectrum,lambda); % Interpolate 
intgreenspectrum = interp1(Osram_lambda,greenspectrum,lambda);
intredspectrum = interp1(Osram_lambda,redspectrum,lambda);

%%

%% New reflectance calculation for LED integrated response

tic

for i = 1:numel(lambda)
for j = 1:numel(Bulkn)
fprintf("Now running lambda %.0f\n",i)
fprintf("Now running ref index %.0f\n",j)

[Refmat1_lam_L0_n(j,i),Z1(j,i)] = multidiel1([Bulkn(j);n_Si(i,2)],0,lambda(i));        

% Calculates reflectance as a function of thickness L from 100 to 120 nm and refractive index n from 1.33 to 1.35 

end
end

toc

Refmat_lam_L0_n = conj(Refmat1_lam_L0_n).*Refmat1_lam_L0_n; %Refmat_Bulk(n,L) 


%% PLOT INTERPOLATED RGB SPECTRUMS

figure(7)
hold on
plot(lambda,intbluespectrum,'b','LineWidth',2) 
plot(lambda,intgreenspectrum,'g','LineWidth',2)
plot(lambda,intredspectrum,'r','LineWidth',2)
xlabel('lambda (nm)')
ylabel('Relative Intensity')
title("RGB Spectrum Data")
saveas(figure(6),[pwd '/Figures/BulkSim/6RGBSpecData.fig']);

%%

for i = 1:numel(Bulkn) 
refcurve0nm_bulk(:,i) = Refmat_lam_L0_n(i,:)';
end

%%

Ref_spec_red_bulk = zeros(3001,21);
Ref_spec_green_bulk = zeros(3001,21);
Ref_spec_blue_bulk = zeros(3001,21);
%%
for i = 1:21
    for j = 1:3001

Ref_spec_red_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intredspectrum); %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intgreenspectrum); %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intbluespectrum); %Reflectivity spectrum for blue at 0 nm (R_b)
    end
end
%%

Ref_spec_red_bulk = Ref_spec_red_bulk';
Ref_spec_green_bulk = Ref_spec_green_bulk';
Ref_spec_blue_bulk = Ref_spec_blue_bulk';
%%


%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)



figure(8)
hold on

plot(lambda,refcurve0nm_bulk(:,find(Bulkn==1.33)),'k','Linewidth',2); %Reflectivity curve at silicon at n = 1.33
plot(lambda,refcurve0nm_bulk(:,find(Bulkn==1.35)),'k','Linewidth',2); %Reflectivity curve at silicon at n = 1.35
plot(lambda,intredspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,intbluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication
plot(lambda,Ref_spec_red_bulk(find(Bulkn==1.33),:),'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(find(Bulkn==1.33),:),'b','Linewidth',2) %Blue spectrum after multiplication
plot(lambda,Ref_spec_red_bulk(find(Bulkn==1.35),:),'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(find(Bulkn==1.35),:),'b','Linewidth',2) %Blue spectrum after multiplication
xlim([400 700])
ylim([0 1.3])
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon (no oxide)')
legend('Reflectivity Curve for Silicon @ n = 1.33','Reflectivity Curve for Silicon @ n = 1.35','Red Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','northeast')
saveas(figure(7),[pwd '/Figures/RefSim/7SiRef.fig']);

%% CALCULATE AGAIN FROM BEGINNING

Refmat_lam_L0_n = Refmat_lam_L0_n'
%%
I_refred_bulk_L0 = []
I_refgreen_bulk_L0 = []
I_refblue_bulk_L0 = []
%
I_refred_bulk = []; %Reflected Intensity for red
I_refgreen_bulk = []; %Reflected Intensity for green
I_refblue_bulk = []; %Reflected Intensity for blue

for i=1:size(Refmat_lam_L_n,2)
    for j = 1:numel(Bulkn)
    fprintf("Now running %.0f\n",i)
    I_refred_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_red_bulk(j,2:end)')./100;
    I_refgreen_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_green_bulk(j,2:end)')./100;
    I_refblue_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_blue_bulk(j,2:end)')./100;
    I_refred_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_red_bulk(j,2:end)')./100;
    I_refgreen_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_green_bulk(j,2:end)')./100;
    I_refblue_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_blue_bulk(j,2:end)')./100;
    end 
end

%%

% hold on
% plot(Bulkn,I_refred_bulk_L0(1,:))
% plot(lambda(2:end),I_refred_bulk_L0(:,1))


%% Plot total reflected intensity for RGB
figure(9)
hold on
plot(BulkL,I_refblue_bulk(:,:,find(Bulkn==1.33)),'b','LineWidth',2)
plot(BulkL,I_refred_bulk(:,:,find(Bulkn==1.33)),'r','LineWidth',2)
plot(BulkL,I_refblue_bulk(:,:,find(Bulkn==1.34)),'b--','LineWidth',2)
plot(BulkL,I_refred_bulk(:,:,find(Bulkn==1.34)),'r--','LineWidth',2)
plot(BulkL,I_refblue_bulk(:,:,find(Bulkn==1.35)),'b-.','LineWidth',2)
plot(BulkL,I_refred_bulk(:,:,find(Bulkn==1.35)),'r-.','LineWidth',2)
title('Reflected Intensity (I_O_u_t) for different refractive indices')
legend('n = 1.33','n = 1.33','n = 1.34','n = 1.34','n = 1.35','n = 1.35')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT AS A FUNCTION OF REF INDEX WITH DIFFERENT THICKNESSES (L = 100, 105,..., 120)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_refred_bulk = reshape(I_refred_bulk,201,21);
I_refgreen_bulk = reshape(I_refgreen_bulk,201,21);
I_refblue_bulk = reshape(I_refblue_bulk,201,21);

%%

figure(10)
hold on
plot(Bulkn,I_refblue_bulk(find(BulkL == 100),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 105),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 110),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 115),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 120),:),'LineWidth',2)
xlim([1.33 1.35])
legend('L=100','L=105','L=110','L=115','L=120')
title("Reflected Intensity for blue")
%%
figure(11)
hold on
plot(Bulkn,I_refred_bulk(find(BulkL == 100),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 105),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 110),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 115),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 120),:),'LineWidth',2)
xlim([1.33 1.35])
legend('L=100','L=105','L=110','L=115','L=120')
title("Reflected Intensity for red")

%% Get ratio

subtrSiT1bulkR = (I_refred_bulk_L0 - I_refred_bulk)' ;
addSiT1bulkR = (I_refred_bulk_L0 + I_refred_bulk)'; 
subtrSiT1bulkB = (I_refblue_bulk_L0 - I_refblue_bulk)' ;
addSiT1bulkB = (I_refblue_bulk_L0 + I_refblue_bulk)';

ratSiT1bulkR = subtrSiT1bulkR ./ addSiT1bulkR;
ratSiT1bulkB = subtrSiT1bulkB ./ addSiT1bulkB;
%% 

figure(14) 

hold on
plot(Bulkn,ratSiT1bulkR(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkR(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkR(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkR(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkR(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')

%% 

figure(15) 

hold on
plot(Bulkn,ratSiT1bulkB(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkB(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkB(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkB(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratSiT1bulkB(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (B)')

%% Get and display slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slpratbulkR(i) = (ratSiT1bulkR(end,i) - ratSiT1bulkR(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

for i = 1:numel(BulkL)
    slpratbulkB(i) = (ratSiT1bulkB(end,i) - ratSiT1bulkB(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

%%
figure(16)
hold on
plot(BulkL,slpratbulkR,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')
ylim([-1.29 -1.22])

figure(17)
hold on
plot(BulkL,slpratbulkB,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (B)')
ylim([-0.85 -0.35])

%% COMPARISON

figure(18)
hold on
plot(BulkL,slpratbulkR,'LineWidth',2)
plot(BulkL,slpratR,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R) difference between SW and LED int')
ylim([-1.29 -1.22])
legend('LED Int','Single Wavelength')

figure(19)
hold on
plot(BulkL,slpratbulkB,'LineWidth',2)
plot(BulkL,slpratB,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (B) difference between SW and LED int')
ylim([-0.85 -0.35])
legend('LED Int','Single Wavelength')
