clc;clear;close all;addpath(genpath(pwd))

load("SimOutputs/RefSimtoBulkData.mat") % Loads variables needed for simulations

%% Get R_Si and R_T1 for air

R_Si_air = Refmat_air(find(lambda==cw_r),find(L==0)); % Reflectance at Si
R_T1_air = Refmat_air(find(lambda==cw_r),find(L==100):find(L==120)); % Reflectance at T1 = 100 nm to 120 nm oxide

%% Get ratio R_Si - R_T1 / R_Si + R_T1

subtrair = R_Si_air - R_T1_air; % Subtract T1 oxide reflectance from Silicon reflectance
addair = R_Si_air + R_T1_air; % Add T1 oxide reflectance and Silicon reflectance
rat = subtrair ./ addair; % Ratio 

%% Display Ratio as a function of thickness

figure(1)
plot(linspace(100,120,201),rat)
xlabel('Thickness T1 (nm)'); ylabel('Ratio')
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1) (FOR AIR)')
saveas(figure(1),[pwd '/Figures/BulkSim/1RatAir.fig']);
%% Define n and L vars

Z1_bulk = [];
Refmat1_bulk = [];

Bulkn = linspace(1.33,1.35,21); % Refractive index from n = 1.33 to 1.35
BulkL = linspace(100,120,201); % Thickness from 100 nm to 120 nm

%% Get reflectance vals from n = 1.33 to 1.35 and silicon for red single wavelength

%% CHANGE N_SIO2 N_SI

for i = 1
    for j = 1
        for k = 1:numel(Bulkn)
        [Refmat_lamR_L0_n_1(i,j,k),Z_lamR_L0_n_1(i,j,k)] = multidiel1([Bulkn(k);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
        end
    end
end
%% Calculates reflectance at Silicon (L = 0) as a function of refractive index (n from 1.33 to 1.35)

Refmat_lamR_L0_n = conj(Refmat_lamR_L0_n_1).*Refmat_lamR_L0_n_1;

%%

for i = 1:numel(Bulkn)
[Refmat1Si_bulk(i),Z1(i)] = multidiel1([Bulkn(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
end

% Calculates reflectance at Silicon (L = 0) as a function of refractive index (n from 1.33 to 1.35)

RefmatSi_bulk = conj(Refmat1Si_bulk).*Refmat1Si_bulk;

%%
Refmat_lam_L_n_1 = []
Z1_lam_L_n_1 = []

tic

for i = 1:numel(lambda)
    for j = 1:numel(BulkL)
        for k = 1:numel(Bulkn)
            fprintf("Now running lambda %.0f thickness %.0f ref index %.0f\n",i,j,k)
            Refmat_lam_L_n_1(i,j,k) = multidiel1([Bulkn(k);n_SiO2(i,2);n_Si(i,2)],BulkL(j).*n_SiO2(i,2),lambda(i));
        end
    end
end

toc


%%
Refmat_lam_L_n = conj(Refmat_lam_L_n_1).*Refmat_lam_L_n_1;

%% Get reflectance vals from n = 1.33 to 1.35 and thickness 100 120 nm for red led

tic

for i = 1:numel(Bulkn)
for j = 1:numel(BulkL)
fprintf("Now running ref index %.0f\n",i)
fprintf("Now running thickness %.0f\n",j)

[Refmat1_bulk(i,j),Z1(i,j)] = multidiel1([Bulkn(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],BulkL(j).*n_SiO2(find(lambda==cw_r),2),cw_r);        

% Calculates reflectance as a function of thickness L from 100 to 120 nm and refractive index n from 1.33 to 1.35 

end
end

toc

%%

Refmat_bulk = conj(Refmat1_bulk).*Refmat1_bulk; %Refmat_Bulk(n,L) 

%% Display reflectance for 5 different thickness from 100 to 120 nm for single wavelength lambda = 633 nm


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
saveas(figure(2),[pwd '/Figures/BulkSim/2Refvs5L.fig']);
%%
% figure(3)
% hold on
% plot(Bulkn,Refmat_bulk(:,find(BulkL==100)),'LineWidth',2)
% plot(Bulkn,Refmat_bulk(:,find(BulkL==105)),'LineWidth',2)
% plot(Bulkn,Refmat_bulk(:,find(BulkL==110)),'LineWidth',2)
% plot(Bulkn,Refmat_bulk(:,find(BulkL==115)),'LineWidth',2)
% plot(Bulkn,Refmat_bulk(:,find(BulkL==120)),'LineWidth',2)
% 
% legend('L=100','L=105','L=110','L=115','L=120')
% xlabel('ref index');ylabel('Reflectance')
% xlim([1.33 1.35])
% title('Reflectance')
% saveas(figure(2),[pwd '/Figures/BulkSim/2Refvs5L.fig']);

%% Get T0-T1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET THESE ACCORDING TO CHANGED VARIABLE NAMES AND FORMAT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Refmat_lamR_L0_n = reshape(Refmat_lamR_L0_n,[21 1]); % Reflectance at silicon for n from 1.33 to 1.35
%%
Refmat_lamR_L_n = reshape(Refmat_lam_L_n(find(lambda == cw_r),:,:),[201 21])';
%%
subtrSiT1 = Refmat_lamR_L0_n - Refmat_lamR_L_n;
addSiT1 = Refmat_lamR_L0_n + Refmat_lamR_L_n;
ratSiT1 = subtrSiT1 ./ addSiT1;
%%

%%
% RefmatSi_bulk = RefmatSi_bulk';
% subbulk = RefmatSi_bulk - Refmat_bulk; % Subtract T1 oxide reflectance from Silicon reflectance
% addbulk = RefmatSi_bulk + Refmat_bulk; % Add T1 oxide reflectance and Silicon reflectance
% ratbulk = subbulk ./ addbulk; % Ratio
%% 



%% Display ratio R_Si - R_T1 / R_Si + R_T1 as a function of n from 1.33 to 1.35

% figure(30) 
% 
% hold on
% plot(Bulkn,ratbulk(:,find(BulkL==100)),'LineWidth',2)
% plot(Bulkn,ratbulk(:,find(BulkL==105)),'LineWidth',2)
% plot(Bulkn,ratbulk(:,find(BulkL==110)),'LineWidth',2)
% plot(Bulkn,ratbulk(:,find(BulkL==115)),'LineWidth',2)
% plot(Bulkn,ratbulk(:,find(BulkL==120)),'LineWidth',2)
% legend('L=100','L=105','L=110','L=115','L=120')
% xlabel('ref index');ylabel('Ratio')
% xlim([1.33 1.35])
% title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
% saveas(figure(3),[pwd '/Figures/BulkSim/3RatBulk.fig']);

%%

figure(3) 

hold on
plot(Bulkn,ratSiT1(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratSiT1(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratSiT1(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratSiT1(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratSiT1(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(3),[pwd '/Figures/BulkSim/3RatBulk.fig']);


%% Get and display slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slprat(i) = (ratSiT1(end,i) - ratSiT1(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end



%% Display slope of ratio as a function of thickness

figure(4)
hold on
plot(BulkL,slprat,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(4),[pwd '/Figures/BulkSim/4SlpRat.fig']);

%% Code to test n in a broader range (1 to 3)

%% BURAYA KADAR OK: BURADAN SONRA CHANGE VAR NAMES

ntest = linspace(1,3,201); 


for i = 1:numel(ntest)
[RefmatSitest_bulk(i),Z1(i)] = multidiel1([ntest(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
end

% Calculates reflectance at Silicon (L = 0) as a function of refractive index (n from 1.33 to 1.35)

RefmatSitest_bulk = conj(RefmatSitest_bulk).*RefmatSitest_bulk;

%% Get reflectance vals from n = 1.33 to 1.35 and thickness 100 120 nm for red led

tic

for i = 1:numel(ntest)
for j = 1:numel(BulkL)
fprintf("Now running ref index %.0f\n",i)
fprintf("Now running thickness %.0f\n",j)

[Refmat1test_bulk(i,j),Z1(i,j)] = multidiel1([ntest(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],BulkL(j).*n_SiO2(find(lambda==cw_r),2),cw_r);        

% Calculates reflectance as a function of thickness L from 100 to 120 nm and refractive index n from 1.33 to 1.35 

end
end

toc

Refmattest_bulk = conj(Refmat1test_bulk).*Refmat1test_bulk; %Refmat_Bulk(n,L) 

RefmatSitest_bulk = RefmatSitest_bulk';
subbulktest = RefmatSitest_bulk - Refmattest_bulk; % Subtract T1 oxide reflectance from Silicon reflectance
addbulktest = RefmatSitest_bulk + Refmattest_bulk; % Add T1 oxide reflectance and Silicon reflectance
ratbulktest = subbulktest ./ addbulktest; % Ratio
%% 

figure(5)
hold on
plot(ntest,ratbulktest(:,find(BulkL==100)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==105)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==110)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==115)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(5),[pwd '/Figures/BulkSim/5RatBulktestclose.fig']);

%%

figure(6)
hold on
plot(ntest,ratbulktest(:,find(BulkL==100)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==105)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==110)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==115)),'LineWidth',2)
plot(ntest,ratbulktest(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(5),[pwd '/Figures/BulkSim/6RatBulktest.fig']);

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

for i = 1:numel(Bulkn)
for j = 1:numel(lambda)
fprintf("Now running ref index %.0f\n",i)
fprintf("Now running lambda %.0f\n",j)

[Refmat1_bulk_LEDint(i,j),Z1(i,j)] = multidiel1([Bulkn(i);n_Si(j,2)],0,lambda(j));        

% Calculates reflectance as a function of thickness L from 100 to 120 nm and refractive index n from 1.33 to 1.35 

end
end

toc

Refmat_bulk_LEDint = conj(Refmat1_bulk_LEDint).*Refmat1_bulk_LEDint; %Refmat_Bulk(n,L) 


%% PLOT INTERPOLATED RGB SPECTRUMS

figure(10)
hold on
plot(lambda,intbluespectrum,'b','LineWidth',2) 
plot(lambda,intgreenspectrum,'g','LineWidth',2)
plot(lambda,intredspectrum,'r','LineWidth',2)
xlabel('lambda (nm)')
ylabel('Relative Intensity')
title("RGB Spectrum Data")
saveas(figure(10),[pwd '/Figures/RefSim/10RGBSpecData.fig']);

%%
for i = 1:size(Refmat_bulk_LEDint,1)
Ref_spec_red_bulk(i,:) = Refmat_bulk_LEDint(i,:).*intredspectrum; %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green_bulk(i,:) = Refmat_bulk_LEDint(i,:).*intgreenspectrum; %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue_bulk(i,:) = Refmat_bulk_LEDint(i,:).*intbluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end
%%
refcurve_n133r = Refmat_bulk_LEDint(find(Bulkn==1.33),:);
refcurve_n135r = Refmat_bulk_LEDint(find(Bulkn==1.35),:);

%% Get a matrix with lambda,L and n
% 
% Refmatwn1 = []
% Zwn1 = []
% 
% [Refmatwn(i,j,k),Z1(i,j,)] = multidiel1([Bulkn(k);n_Air(j,2);n_Si(j,2)],0,lambda(j));



%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)
figure(11)
hold on

plot(lambda,refcurve_n133r,'k','Linewidth',2);
plot(lambda,refcurve_n135r,'m','Linewidth',2);
plot(lambda,intredspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,intbluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication

plot(lambda,Ref_spec_red_bulk(1,:),'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(1,:),'b','Linewidth',2) %Blue spectrum after multiplication
plot(lambda,Ref_spec_red_bulk(21,:),'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_blue_bulk(21,:),'b','Linewidth',2) %Blue spectrum after multiplication
xlim([400 700])
ylim([0 1.3])
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon (no oxide)')
legend('Reflectivity Curve for Silicon @ n = 1.33','Reflectivity Curve for Silicon @ n = 1.35','Red Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','northeast')
saveas(figure(11),[pwd '/Figures/RefSim/11SiRef.fig']);

%% Calculate total reflected intensity by multiplying reflectance with reflected intensity (I_Out)

I_refred_bulk = []; %Reflected Intensity for red
I_refgreen_bulk = []; %Reflected Intensity for green
I_refblue_bulk = []; %Reflected Intensity for blue
for i=1:size(Refmat_bulk_LEDint,2)
   fprintf("Now running %.0f\n",i)
   I_refred_bulk(:,i) = sum(Refmat_bulk_LEDint(2:end,i).*Ref_spec_red_bulk(2:end,i));      
   I_refgreen_bulk(:,i) = sum(Refmat_bulk_LEDint(2:end,i).*Ref_spec_green_bulk(2:end,i));     
   I_refblue_bulk(:,i) = sum(Refmat_bulk_LEDint(2:end,i).*Ref_spec_blue_bulk(2:end,i));       
end
%% Plot total reflected intensity for RGB
figure(12)
hold on
plot(L,I_refblue/100,'b','LineWidth',2)
plot(L,I_refgreen/100,'g','LineWidth',2)
plot(L,I_refred/100,'r','LineWidth',2)
title('Reflected Intensity (I_O_u_t) for RGB from 0 nm to 300 nm')
legend('Blue','Green','Red')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);
saveas(figure(12),[pwd '/Figures/RefSim/12RefISiO2.fig']);





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





%%









