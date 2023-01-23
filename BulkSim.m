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

%% Get reflectance vals from n = 1.33 to 1.35 and silicon for red led

for i = 1:numel(Bulkn)
[Refmat1Si_bulk(i),Z1(i)] = multidiel1([Bulkn(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r); 
end

% Calculates reflectance at Silicon (L = 0) as a function of refractive index (n from 1.33 to 1.35)

RefmatSi_bulk = conj(Refmat1Si_bulk).*Refmat1Si_bulk;

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

%% Display reflectance for 5 different thickness from 100 to 120 nm

figure(2)
hold on
plot(Bulkn,Refmat_bulk(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,Refmat_bulk(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,Refmat_bulk(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,Refmat_bulk(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,Refmat_bulk(:,find(BulkL==120)),'LineWidth',2)

legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Reflectance')
xlim([1.33 1.35])
title('Reflectance')
saveas(figure(2),[pwd '/Figures/BulkSim/2Refvs5L.fig']);

%% Get T0-T1

RefmatSi_bulk = RefmatSi_bulk';
subbulk = RefmatSi_bulk - Refmat_bulk; % Subtract T1 oxide reflectance from Silicon reflectance
addbulk = RefmatSi_bulk + Refmat_bulk; % Add T1 oxide reflectance and Silicon reflectance
ratbulk = subbulk ./ addbulk; % Ratio
%% 



%% Display ratio R_Si - R_T1 / R_Si + R_T1 as a function of n from 1.33 to 1.35

figure(3) 

hold on
plot(Bulkn,ratbulk(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn,ratbulk(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn,ratbulk(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn,ratbulk(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn,ratbulk(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index');ylabel('Ratio')
xlim([1.33 1.35])
title ('(R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(3),[pwd '/Figures/BulkSim/3RatBulk.fig']);

%% Get and display slope of ratio for L from 100 to 120 nm

for i = 1:numel(BulkL)
    slprat(i) = (ratbulk(end,i) - ratbulk(1,i)) / (1.35-1.33); % Calculates slope of ratio R_Si - R_T1 / R_Si + R_T1
end



%% Display slope of ratio as a function of thickness

figure(4)
hold on
plot(BulkL,slprat,'LineWidth',2)
xlabel('L (nm)'); ylabel('Slope of Ratio');
xlim([100 120])
title ('Slope of ratio (R_S_i - R_T_1) / (R_S_i + R_T_1)')
saveas(figure(4),[pwd '/Figures/BulkSim/4SlpRat.fig']);












%% Extra Code (TEST n from 1 -> 1,33)

%ntest = linspace(1:1.35)
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









