clc;clear;close all;addpath(genpath(pwd))
%%
load("SimOutputs/RefSimData.mat")

%% For air

%% Get R_Si and R_T1

R_Si_air = Refmat_air(find(lambda==cw_r),find(L==0));
R_T1_air = Refmat_air(find(lambda==cw_r),find(L==100):find(L==120));
%%
subtr = R_Si_air - R_T1_air;
add = R_Si_air + R_T1_air;
rat = subtr ./ add

%% Display Ratio

figure(1)
plot(linspace(100,120,201),rat)
xlabel('Thickness T1 (nm'); ylabel('Ratio')
title ('R_S_i - R_T_1 / R_S_i + R_T_1 (FOR AIR)')

%% Get n from 1.33 to 1.35



%%

% for i = 1:21
% Bulkn(:,1,i) = lambda;
% Bulkn(:,2,i) = 1.33+0.001*(i-1);
% end


%%

%% GET REFLECTANCE EQUATION FOR SINGLE COLOR

Z1_bulk = [];
Refmat1_bulk = [];
%%
Bulkn = linspace(1.33,1.35,21);
BulkL = linspace(100,120,201);

%% Get reflectance vals from n = 1.33 to 1.35 and silicon for red led

for i = 1:numel(Bulkn)
[Refmat1Si_bulk(i),Z1(i)] = multidiel1([Bulkn(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],0,cw_r)         ;
end

RefmatSi_bulk = conj(Refmat1Si_bulk).*Refmat1Si_bulk;

%% Get reflectance vals from n = 1.33 to 1.35 and thickness 100 120 nm for red led

tic
for i = 1:numel(Bulkn)
for j = 1:numel(BulkL)
fprintf("Now running ref index %.0f\n",i)
fprintf("Now running thickness %.0f\n",j)

[Refmat1_bulk(i,j),Z1(i,j)] = multidiel1([Bulkn(i);n_SiO2(find(lambda==cw_r),2);n_Si(find(lambda==cw_r),2)],BulkL(j).*n_SiO2(find(lambda==cw_r),2),cw_r);       

% tic
% for i = 1:numel(lambda)
% for j = 1:numel(BulkL)
% for k = 1:size(Bulkn,3)
% fprintf("Now running lambda %.0f\n",i)
% fprintf("Now running thickness %.0f\n",j)
% fprintf("Now running ref index %.0f\n",k)
% 
% [Refmat1_bulk(i,j,k),Z1(i,j,k)] = multidiel1([Bulkn(:,2,k);n_SiO2(i,2);n_Si(i,2)],BulkL(j).*n_SiO2(i,2),lambda(i));

end
end

toc
%%

Refmat_bulk = conj(Refmat1_bulk).*Refmat1_bulk; %Refmat_Bulk(n,L) 
%%
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

%% Get T0-T1

RefmatSi_bulk = RefmatSi_bulk';
subbulk = RefmatSi_bulk - Refmat_bulk;
addbulk = RefmatSi_bulk + Refmat_bulk;
ratbulk = subbulk ./ addbulk
%% 



%%

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
title ('R_S_i - R_T_1 / R_S_i + R_T_1')

%%

slprat = diff(ratbulk)./0.001;


%%

figure(4) 

hold on
plot(Bulkn(2:end),slprat(:,find(BulkL==100)),'LineWidth',2)
plot(Bulkn(2:end),slprat(:,find(BulkL==105)),'LineWidth',2)
plot(Bulkn(2:end),slprat(:,find(BulkL==110)),'LineWidth',2)
plot(Bulkn(2:end),slprat(:,find(BulkL==115)),'LineWidth',2)
plot(Bulkn(2:end),slprat(:,find(BulkL==120)),'LineWidth',2)
legend('L=100','L=105','L=110','L=115','L=120')
xlabel('ref index'); ylabel('Slope of Reflectance');
xlim([1.33 1.35])
title ('Slope')



%%

figure(5)
hold on
plot(BulkL, slprat (1, : ))
plot(BulkL, slprat (5, : ))
plot(BulkL, slprat (10, :) )
plot(BulkL, slprat (15, : ) )
plot(BulkL, slprat (20, : ))
%%
% figure(6)
% hold on
% plot(Bulkn(2: end) , slprat(:, find (BulkL==105)) , 'LineWidth',2)
% plot(Bulkn (2: end) , slprat(:,find (BulkL==110)), 'LineWidth',2)
% plot(Bulkn (2: end) , slprat(:,find (BulkL==115)), 'LineWidth',2)
% plot(Bulkn(2: end) , slprat(:, find (BulkL==120)), 'LineWidth', 2)





%%









