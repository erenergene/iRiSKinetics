%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARES TWO DIFFERENT LED INTEGRATION METHOD %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Get reflectance curve at Silicon

for i = 1:numel(Bulkn) 
refcurve0nm_bulk(:,i) = Refmat_lam_L0_n(i,:)';
end

%% Calculate reflectivity spectrum for RGB

Ref_spec_red_bulk = zeros(3001,21);
Ref_spec_green_bulk = zeros(3001,21);
Ref_spec_blue_bulk = zeros(3001,21);

for i = 1:21
    for j = 1:3001

    Ref_spec_red_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intredspectrum); %Reflectivity spectrum for red at 0 nm (R_r)
    Ref_spec_green_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intgreenspectrum); %Reflectivity spectrum for green at 0 nm (R_g)
    Ref_spec_blue_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intbluespectrum); %Reflectivity spectrum for blue at 0 nm (R_b)
    
    end
end

Ref_spec_red_bulk = Ref_spec_red_bulk';
Ref_spec_green_bulk = Ref_spec_green_bulk';
Ref_spec_blue_bulk = Ref_spec_blue_bulk';

%% FIRST METHOD (BulkSim.m)

%% Calculate reflectivity spectrum for RGB

Ref_spec_red_bulk = zeros(3001,21);
Ref_spec_green_bulk = zeros(3001,21);
Ref_spec_blue_bulk = zeros(3001,21);

for i = 1:21
    for j = 1:3001

    Ref_spec_red_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intredspectrum); %Reflectivity spectrum for red at 0 nm (R_r)
    Ref_spec_green_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intgreenspectrum); %Reflectivity spectrum for green at 0 nm (R_g)
    Ref_spec_blue_bulk(:,i) = (refcurve0nm_bulk(j,i)'.*intbluespectrum); %Reflectivity spectrum for blue at 0 nm (R_b)
    
    end
end

Ref_spec_red_bulk = Ref_spec_red_bulk';
Ref_spec_green_bulk = Ref_spec_green_bulk';
Ref_spec_blue_bulk = Ref_spec_blue_bulk';

%% SECOND METHOD (BulkSimPast.m)

for i = 1:size(Refmat_lam_L0_n,1)
Ref_spec_red_bulk(i,:) = Refmat_lam_L0_n(i,:).*intredspectrum; %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green_bulk(i,:) = Refmat_lam_L0_n(i,:).*intgreenspectrum; %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue_bulk(i,:) = Refmat_lam_L0_n(i,:).*intbluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end
%%
refcurve_n133r = Refmat_lam_L0_n(find(Bulkn==1.33),:);
refcurve_n135r = Refmat_lam_L0_n(find(Bulkn==1.35),:);



























