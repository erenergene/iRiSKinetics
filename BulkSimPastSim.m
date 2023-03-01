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
legend('red','green','blue')
saveas(figure(1),[pwd '/Figures/BulkSim/1RatAir.fig']);
saveas(figure(1),[pwd '/Figures/BulkSimJpg/1RatAir.jpg']);


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
title('Reflectance for Red SW')
saveas(figure(2),[pwd '/Figures/BulkSim/2RefRvs5L.fig']);
saveas(figure(2),[pwd '/Figures/BulkSimJpg/2RefRvs5L.jpg']);

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
title('Reflectance for Blue SW')
saveas(figure(3),[pwd '/Figures/BulkSim/3RefBvs5L.fig']);
saveas(figure(3),[pwd '/Figures/BulkSimJpg/3RefBvs5L.jpg']);

%% Get T0 and T1 - %% Reconstruct data for easier calculation

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
saveas(figure(4),[pwd '/Figures/BulkSim/4RatswR5Lvsn.fig']);
saveas(figure(4),[pwd '/Figures/BulkSimJpg/4RatswR5Lvsn.jpg']);
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
saveas(figure(5),[pwd '/Figures/BulkSim/5RatswB5Lvsn.fig']);
saveas(figure(5),[pwd '/Figures/BulkSimJpg/5RatswB5Lvsn.jpg']);

%% Get and display slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slpratR(i) = (ratSiT1R(end,i) - ratSiT1R(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

for i = 1:numel(BulkL)
    slpratB(i) = (ratSiT1B(end,i) - ratSiT1B(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end

%% Display slope of ratio as a function of thickness for red single wavelength

figure(6)
hold on
plot(BulkL,slpratR,'r','LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')
saveas(figure(6),[pwd '/Figures/BulkSim/6SlpRatswRvsL.fig']);
saveas(figure(6),[pwd '/Figures/BulkSimJpg/6SlpRatswRvsL.jpg']);


%% Display slope of ratio as a function of thickness for blue single wavelength

figure(7)
hold on
plot(BulkL,slpratB,'b','LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (B)')
saveas(figure(7),[pwd '/Figures/BulkSim/7SlpRatswBvsL.fig']);
saveas(figure(7),[pwd '/Figures/BulkSimJpg/7SlpRatswBvsL.jpg']);

%% Code to visualize n in a broader range (n = 1 to 3)

nbroad = linspace(1,3,201); %% Define broad n

%% Get reflectance vals: lambda = 400 - 700 nm, thickness = 100 - 120 nm, ref index = 1 - 3


for i = 1:numel(nbroad)
[Refmat1_lamR_L0_nbroad(i),Z1(i)] = multidiel1([nbroad(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
end

Refmat_lamR_L0_nbroad = conj(Refmat1_lamR_L0_nbroad).*Refmat1_lamR_L0_nbroad;

%% Get reflectance vals: lambda = 400 - 700 nm, thickness = 0 nm, ref index = 1 - 3

tic

for i = 1:numel(nbroad)
for j = 1:numel(BulkL)
[Refmat1_lamR_L_nbroad(i,j),Z1(i,j)] = multidiel1([nbroad(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],BulkL(j).*n_SiO2(find(lambda==cw_r),2),cw_r);        
end
end

toc

Refmat_lamR_L_nbroad = conj(Refmat1_lamR_L_nbroad).*Refmat1_lamR_L_nbroad;

Refmat_lamR_L_nbroad = Refmat_lamR_L_nbroad';
subbulknbroad = Refmat_lamR_L0_nbroad - Refmat_lamR_L_nbroad; % Subtract T1 oxide reflectance from Silicon reflectance
addbulknbroad = Refmat_lamR_L0_nbroad + Refmat_lamR_L_nbroad; % Add T1 oxide reflectance and Silicon reflectance
ratbulknbroad = (subbulknbroad ./ addbulknbroad)'; % Ratio

%% Plot broad n vs ratio

figure(8)
hold on
plot(nbroad,ratbulknbroad(:,find(BulkL==100)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==105)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==110)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==115)),'LineWidth',2)
plot(nbroad,ratbulknbroad(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(8),[pwd '/Figures/BulkSim/8Ratbroadntest.fig']);
saveas(figure(8),[pwd '/Figures/BulkSimJpg/8Ratbroadntest.jpg']);

%% LED Integrated Response

load("Osram_Spec_Data_edited.mat") % loads LED data
bluespectrum = Spec_DentalBlue_460nm; clear("Spec_DentalBlue_460nm")
greenspectrum = Spec_Green_517nm; clear("Spec_Green_517nm")
redspectrum = Spec_Red_633nm; clear("Spec_Red_633nm")

intbluespectrum = interp1(Osram_lambda,bluespectrum,lambda); % Interpolate 
intgreenspectrum = interp1(Osram_lambda,greenspectrum,lambda);
intredspectrum = interp1(Osram_lambda,redspectrum,lambda);

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

figure(9)
hold on
plot(lambda,intbluespectrum,'b','LineWidth',2) 
plot(lambda,intgreenspectrum,'g','LineWidth',2)
plot(lambda,intredspectrum,'r','LineWidth',2)
xlabel('lambda (nm)')
ylabel('Relative Intensity')
title("RGB Spectrum Data")
saveas(figure(9),[pwd '/Figures/BulkSim/9RGBSpecData.fig']);
saveas(figure(9),[pwd '/Figures/BulkSimJpg/9RGBSpecData.jpg']);

%%
Ref_spec_red_n1 = Refmat_air(:,find(L==0))'.*intredspectrum;
Ref_spec_green_n1 = Refmat_air(:,find(L==0))'.*intgreenspectrum;
Ref_spec_blue_n1 = Refmat_air(:,find(L==0))'.*intbluespectrum;
%%
for i = 1:size(Refmat_lam_L0_n,1)
Ref_spec_red_bulk(i,:) = Refmat_lam_L0_n(i,:).*intredspectrum; %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green_bulk(i,:) = Refmat_lam_L0_n(i,:).*intgreenspectrum; %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue_bulk(i,:) = Refmat_lam_L0_n(i,:).*intbluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end
%%
refcurve_n1r = refcurve0nm
refcurve_n133r = Refmat_lam_L0_n(find(Bulkn==1.33),:);
refcurve_n135r = Refmat_lam_L0_n(find(Bulkn==1.35),:);

%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)
figure(10)
hold on

plot(lambda,refcurve_n1r,'c','Linewidth',2);
plot(lambda,refcurve_n133r,'k','Linewidth',2);
plot(lambda,refcurve_n135r,'m','Linewidth',2);
plot(lambda,intredspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,intgreenspectrum,'g','Linewidth',2) %Blue spectrum before multiplication
plot(lambda,intbluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication

plot(lambda,Ref_spec_red_n1,'c','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green_n1,'c','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_n1,'c','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_red_bulk(1,:),'k','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green_bulk(1,:),'k','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(1,:),'k','Linewidth',2) %Blue spectrum after multiplication
plot(lambda,Ref_spec_red_bulk(21,:),'m','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green_bulk(21,:),'m','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(21,:),'m','Linewidth',2) %Blue spectrum after multiplication
xlim([400 700])
ylim([0 1.3])
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon (no oxide)')
legend('Reflectivity Curve for Silicon @ n = 1.33','Reflectivity Curve for Silicon @ n = 1.35','Red Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','northeast')
saveas(figure(10),[pwd '/Figures/BulkSimPast/7SiRef.fig']);
%%

%% Calculate reflected intensity for L = 0 and L = 100 - 120 and n = 1.33 -1.35

Refmat_lam_L0_n = Refmat_lam_L0_n'

I_refred_bulk_L0 = []
I_refgreen_bulk_L0 = []
I_refblue_bulk_L0 = []

I_refred_bulk = []; %Reflected Intensity for red
I_refgreen_bulk = []; %Reflected Intensity for green
I_refblue_bulk = []; %Reflected Intensity for blue

for i=1:size(Refmat_lam_L_n,2)
    for j = 1:numel(Bulkn)
    fprintf("Now running %.0f\n",i)
    I_refred_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_red_bulk(j,2:end)')./100; %Red - L = 0
    I_refgreen_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_green_bulk(j,2:end)')./100; %Green - L = 0
    I_refblue_bulk_L0(:,j) = sum(Refmat_lam_L0_n(2:end,j).*Ref_spec_blue_bulk(j,2:end)')./100; %Blue - L = 0
    I_refred_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_red_bulk(j,2:end)')./100; %Red - L = 100-120
    I_refgreen_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_green_bulk(j,2:end)')./100; %Green - L = 100-120
    I_refblue_bulk(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_blue_bulk(j,2:end)')./100; %Blue - L = 100-120
    end 
end

%% Plot total reflected intensity for R and B for n = 1.33, 1.34, and 1.35 as a function of thickness

figure(11)
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
saveas(figure(11),[pwd '/Figures/BulkSim/11RefIntRB3nvsL.fig']);
saveas(figure(11),[pwd '/Figures/BulkSimJpg/11RefIntRB3nvsL.jpg']);

%% Reconstruct data for easier plot

I_refred_bulk = reshape(I_refred_bulk,201,21);
I_refgreen_bulk = reshape(I_refgreen_bulk,201,21);
I_refblue_bulk = reshape(I_refblue_bulk,201,21);

%% Plot reflected intensity for 5 thicknesses as a function of n for red

figure(12)
hold on
plot(Bulkn,I_refred_bulk(find(BulkL == 100),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 105),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 110),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 115),:),'LineWidth',2)
plot(Bulkn,I_refred_bulk(find(BulkL == 120),:),'LineWidth',2)
xlim([1.33 1.35])
legend('L=100','L=105','L=110','L=115','L=120')
title("Reflected Intensity for red")
saveas(figure(12),[pwd '/Figures/BulkSim/12RefIntR5Lvsn.fig']);
saveas(figure(12),[pwd '/Figures/BulkSimJpg/12RefIntR5Lvsn.jpg']);

%% Plot reflected intensity for 5 thicknesses as a function of n for blue

figure(13)
hold on
plot(Bulkn,I_refblue_bulk(find(BulkL == 100),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 105),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 110),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 115),:),'LineWidth',2)
plot(Bulkn,I_refblue_bulk(find(BulkL == 120),:),'LineWidth',2)
xlim([1.33 1.35])
legend('L=100','L=105','L=110','L=115','L=120')
title("Reflected Intensity for blue")
saveas(figure(13),[pwd '/Figures/BulkSim/13RefIntB5Lvsn.fig']);
saveas(figure(13),[pwd '/Figures/BulkSimJpg/13RefIntB5Lvsn.jpg']);



%% Get ratio for red LED int

subtrSiT1bulkR = (I_refred_bulk_L0 - I_refred_bulk)' ;
addSiT1bulkR = (I_refred_bulk_L0 + I_refred_bulk)'; 
ratSiT1bulkR = subtrSiT1bulkR ./ addSiT1bulkR;

%% Get ratio for blue LED int

subtrSiT1bulkB = (I_refblue_bulk_L0 - I_refblue_bulk)' ;
addSiT1bulkB = (I_refblue_bulk_L0 + I_refblue_bulk)';
ratSiT1bulkB = subtrSiT1bulkB ./ addSiT1bulkB;

%% Plot ratio for red LED int

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
saveas(figure(14),[pwd '/Figures/BulkSim/14RatIntR5Lvsn.fig']);
saveas(figure(14),[pwd '/Figures/BulkSimJpg/14RatIntR5Lvsn.jpg']);

%% Plot ratio for blue LED int 

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
saveas(figure(15),[pwd '/Figures/BulkSim/15RatIntB5Lvsn.fig']);
saveas(figure(15),[pwd '/Figures/BulkSimJpg/15RatIntB5Lvsn.jpg']);

%% Get slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slpratbulkR(i) = (ratSiT1bulkR(end,i) - ratSiT1bulkR(1,i)) / (1.35-1.33); % Calculates slope of ratio (R_Si - R_T1) / (R_Si + R_T1)
end

for i = 1:numel(BulkL)
    slpratbulkB(i) = (ratSiT1bulkB(end,i) - ratSiT1bulkB(1,i)) / (1.35-1.33); % Calculates slope of ratio (R_Si - R_T1) / (R_Si + R_T1)
end

%% Display slope of ratio as a function of thickness for red LED integrated response
figure(16)
hold on
plot(BulkL,slpratbulkR,'r','LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R)')
saveas(figure(16),[pwd '/Figures/BulkSim/16SlpRatIntRBvsL.fig']);
saveas(figure(16),[pwd '/Figures/BulkSimJpg/16SlpRatIntRBvsL.jpg']);

%% Display slope of ratio as a function of thickness for red LED integrated response

figure(17)
hold on
plot(BulkL,slpratbulkB,'b','LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (B)')
saveas(figure(17),[pwd '/Figures/BulkSim/17SlpRatIntBvsL.fig']);
saveas(figure(17),[pwd '/Figures/BulkSimJpg/17SlpRatIntBvsL.jpg']);

%% Graph to compare single wavelength and LED integrated response for red

figure(18)
hold on
plot(BulkL,slpratbulkR,'LineWidth',2)
plot(BulkL,slpratR,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (R) difference between SW and LED int')
ylim([-1.29 -1.22])
legend('LED Int','Single Wavelength')
saveas(figure(18),[pwd '/Figures/BulkSim/18SlpRatInt&swRvsL.fig']);
saveas(figure(18),[pwd '/Figures/BulkSimJpg/18SlpRatInt&swRvsL.jpg']);
%% Graph to compare single wavelength and LED integrated response for red

figure(19)
hold on
plot(BulkL,slpratbulkB,'LineWidth',2)
plot(BulkL,slpratB,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1) (B) difference between SW and LED int')
ylim([-0.85 -0.35])
legend('LED Int','Single Wavelength')
saveas(figure(19),[pwd '/Figures/BulkSim/19SlpRatInt&swBvsL.fig']);
saveas(figure(19),[pwd '/Figures/BulkSimJpg/19SlpRatInt&swBvsL.jpg']);


%%
% 
[figures, pathname,]=uigetfile([pwd '/Figures'],'*.jpg','Multiselect','on');
























































%% CALCULATE AGAIN FROM BEGINNING

% I_refred_bulk = []; %Reflected Intensity for red
% I_refgreen_bulk = []; %Reflected Intensity for green
% I_refblue_bulk = []; %Reflected Intensity for blue
% 
% for i=1:size(Refmat_lam_L_n,2)
%     for j = 1:numel(Bulkn)
%     fprintf("Now running %.0f\n",i)
%     I_refred_bulk_red(:,i,j) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_red_bulk(j,2:end)');
%     end 
% end


%% Calculate total reflected intensity by multiplying reflectance with reflected intensity (I_Out)

% I_refred_bulk = []; %Reflected Intensity for red
% I_refgreen_bulk = []; %Reflected Intensity for green
% I_refblue_bulk = []; %Reflected Intensity for blue
% 
% %%
% 
% for i=1:numel(BulkL)
%    fprintf("Now running %.0f\n",i)
%    I_refred_133(:,i) = sum(Refmat_lam_L_n(2:end,i,1).*Ref_spec_red_bulk(1,2:end)');      
%    I_refgreen_133(:,i) = sum(Refmat_lam_L_n(2:end,i,1).*Ref_spec_green_bulk(1,2:end)');         
%    I_refblue_133(:,i) = sum(Refmat_lam_L_n(2:end,i,1).*Ref_spec_blue_bulk(1,2:end)');          
% end
% 
% %%
% 
% for i=1:size(Refmat_lam_L_n,2)
%     for j =1:numel(Bulkn)
%        fprintf("Now running thickness %.0f and ref index %.0f \n",i,j)
%        I_refred_bulk(j,i) = sum(2:end,i,j).*Ref_spec_red_bulk(j,2:end)';
%        I_refgreen_bulk(j,i) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_green_bulk(j,2:end)';
%        I_refblue_bulk(j,i) = sum(Refmat_lam_L_n(2:end,i,j).*Ref_spec_blue_bulk(j,2:end)';
%      end
% end

%% basdan basla




%%
    

%% Plot total reflected intensity for RGB
% figure(12)
% hold on
% plot(L,I_refblue_bulk/100,'b','LineWidth',2)
% plot(L,I_refgreen_bulk/100,'g','LineWidth',2)
% plot(L,I_refred_bulk/100,'r','LineWidth',2)
% title('Reflected Intensity (I_O_u_t) for RGB from 0 nm to 300 nm')
% legend('Blue','Green','Red')
% xlabel('L (\mum)','FontSize',16);
% ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);
% saveas(figure(12),[pwd '/Figures/RefSim/12RefISiO2.fig']);





% figure(4) 
% 
% hold on
% plot(Bulkn(2:end),slprat(:,find(BulkL==100)),'LineWidth',2)
% plot(Bulkn(2:end),slprat(:,find(BulkL==105)),'LineWidth',2)
% plot(Bulkn(2:end),slprat(:,find(BulkL==110)),'LineWidth',2)
% plot(Bulkn(2:end),slprat(:,find(BulkL==115)),'LineWidth',2)
% plot(Bulkn(2:end),slprat(:,find(BulkL==120)),'LineWidth',2)
% legend('L=100','L=105','L=110','L=115','L=120')
% xlabel('ref index'); ylabel('Slope of Reflectance');
% title ('Slope')



%%
% 
% figure(5)
% hold on
% plot(BulkL, slprat (1, : ))
% plot(BulkL, slprat (5, : ))
% plot(BulkL, slprat (10, :) )
% plot(BulkL, slprat (15, : ) )
% plot(BulkL, slprat (20, : ))
%%
% figure(6)
% hold on
% plot(Bulkn(2: end) , slprat(:, find (BulkL==105)) , 'LineWidth',2)
% plot(Bulkn (2: end) , slprat(:,find (BulkL==110)), 'LineWidth',2)
% plot(Bulkn (2: end) , slprat(:,find (BulkL==115)), 'LineWidth',2)
% plot(Bulkn(2: end) , slprat(:, find (BulkL==120)), 'LineWidth', 2)

%FEB1TEST




%%

