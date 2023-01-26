%% WRITE DOWN PSEUDO CODE TOMORROW AT LAB NOTEBOOK

clc; clearvars; close all;
addpath(genpath(pwd))

%% DEFINE WAVELENGTH AND THICKNESS

lambda = 400:0.1:700;                                   %Wavelength: 400-700 nm                          
numval = numel(lambda);
L = linspace(0,300,numval);                             %Thickness: 0-300 nm

cw_b = 460;
cw_g = 517;
cw_r = 633;

%% LOAD REFRACTIVE INDEX VALUES

n_Si_raw = load("Si.nnn");                              %Refractive index of Silicon 
n_SiO2_raw = load("SiO2.nnn");                          %Refractive index of Silicon Oxide

r_Si_raw = load("Si.rrr");                              %Imaginary part of refractive index of Silicon 
r_SiO2_raw = load("SiO2.rrr");                          %Imaginary part of refractive index of Silicon Oxide

k_Si_raw = load("Si.kkk");                              %Extinction Coefficient of Silicon 
k_SiO2_raw = load("SiO2.kkk");                          %Extinction Coefficient of Silicon Oxide

n_Air_raw = load("Borzsonyi.csv");                      %Refractive index of Air
n_Air_raw(:,1) = 1000.* n_Air_raw(:,1);               

n_Water(:,1) = lambda;                                  %Refractive index for solution
n_Water(:,2) = 1.33; 

n_Sol(:,1) = lambda;                                    %Refractive index for solution
n_Sol(:,2) = 1.34; 

%% INTERPOLATE

n_Si = [lambda;interp1(n_Si_raw(:,1),n_Si_raw(:,2),lambda)]';
n_SiO2 = [lambda;interp1(n_SiO2_raw(:,1),n_SiO2_raw(:,2),lambda)]';
n_Air = [lambda;interp1(n_Air_raw(:,1),n_Air_raw(:,2),lambda)]';

%% PLOT INTERPOLATED AND ORIGINAL REFRACTIVE INDICES

figure(1)

subplot(2,2,1) %Refractive Index of Silicon
hold on
plot(n_Si_raw(:,1),n_Si_raw(:,2),"LineWidth",2,'color','b');
plot(lambda,n_Si(:,2),"LineWidth",2,'color','r');
xlabel("Wavelength (nm)")
ylabel("Refractive Index")
title("Si")

subplot(2,2,2) %Refractive index of Silicon Oxide
hold on
plot(n_SiO2_raw(1:58,1),n_SiO2_raw(1:58,2),"LineWidth",2,'color','b');
plot(lambda,n_SiO2(:,2),"LineWidth",2,'color','r');
xlabel("Wavelength (nm)")
ylabel("Refractive Index")
title("SiO2")

subplot(2,2,3) %Refractive index of Air
hold on
plot(n_Air_raw(:,1),n_Air_raw(:,2),"LineWidth",2,'color','b');
plot(lambda,n_Air(:,2),"LineWidth",2,'color','r');
xlabel("Wavelength (nm)")
ylabel("Refractive Index")
title("Air")

subplot(2,2,4) %Interpolated values of refractive index from 400-700nm
hold on
plot(lambda,n_Si(:,2),"LineWidth",2,'color','b');
plot(lambda,n_SiO2(:,2),"LineWidth",2,'color','r');
plot(lambda,n_Air(:,2),"LineWidth",2,'color','g');
xlabel("Wavelength (nm)")
ylabel("Refractive Index")
legend("Si","SiO2","Air")
title('Si & SiO2 for 400-700 nm')
ylim([0 6])

sgtitle('Interpolated vs Raw Refractive Indices')
saveas(figure(1),[pwd '/Figures/RefSim/1RefIndicesIntvsRaw.fig']);


%% CALCULATE REFLECTANCE VALUES

Z1_air = [];
Refmat1_air = [];

tic
for i = 1:numel(lambda)
fprintf("Now running %.0f\n",i)
for j = 1:numel(L)
[Refmat1_air(i,j),Z1(i,j)] = multidiel1([n_Air(i,2);n_SiO2(i,2);n_Si(i,2)],L(j).*n_SiO2(i,2),lambda(i));
[Refmat1_water(i,j),Z1(i,j)] = multidiel1([n_Water(i,2);n_SiO2(i,2);n_Si(i,2)],L(j).*n_SiO2(i,2),lambda(i));
[Refmat1_sol(i,j),Z1(i,j)] = multidiel1([n_Sol(i,2);n_SiO2(i,2);n_Si(i,2)],L(j).*n_SiO2(i,2),lambda(i));
end
end
toc

Refmat_air = conj(Refmat1_air).*Refmat1_air;                
Refmat_water = conj(Refmat1_water).*Refmat1_water;             %Multiply Gamma with conjugate to get rid of imaginary component
Refmat_sol = conj(Refmat1_sol).*Refmat1_sol;

%% DISPLAY REFLECTANCE MATRIX UNDER AIR (n=1)

figure(2)

surf(L,lambda,Refmat_air,'EdgeColor','none');
title('Si/SiO2 Reflectance Under Air')
xlim([0 300])
ylim([400 700])
xlabel('L (\mum)');
ylabel('Lambda (\mum)');
zlabel('Reflectance');
cb = colorbar;
cb.Location = 'eastoutside';
caxis([0 0.5])
saveas(figure(2),[pwd '/Figures/RefSim/2RefMatrixAir.fig']);


%% DISPLAY REFLECTANCE MATRIX UNDER WATER (n=1.33)

figure(3)
surf(L,lambda,Refmat_water,'EdgeColor','none');
title('Si/SiO2 Reflectance Under Water')
xlim([0 300])
ylim([400 700])
xlabel('L (\mum)');
ylabel('Lambda (\mum)');
zlabel('Reflectance');
cb = colorbar;
cb.Location = 'eastoutside';
caxis([0 0.5])
saveas(figure(3),[pwd '/Figures/RefSim/3RefMatrixWater.fig']);


%% THICKNESS VS REFLECTANCE FOR RGB UNDER AIR (n=1)

figure(4)
hold on
plot(L,Refmat_air(find(lambda == cw_b),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L,Refmat_air(find(lambda == cw_g),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L,Refmat_air(find(lambda == cw_r),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity Under Air')
xlabel('L (nm)','FontSize',16);
ylabel('Reflectivity ','FontSize',16');
saveas(figure(4),[pwd '/Figures/RefSim/4LvsRefAirRGB.fig']);


%% THICKNESS VS REFLECTANCE FOR RGB UNDER WATER (n=1.33)

figure(5)
hold on
plot(L,Refmat_water(find(lambda == cw_b),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L,Refmat_water(find(lambda == cw_g),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L,Refmat_water(find(lambda == cw_r),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity Under Water')
xlabel('L (nm)','FontSize',16);
ylabel('Reflectivity ','FontSize',16');
saveas(figure(5),[pwd '/Figures/RefSim/5LvsRefWaterRGB.fig']);


%% THICKNESS VS REFLECTANCE FOR RGB UNDER SOLUTION (n=1.34)

figure(6)
hold on
plot(L,Refmat_sol(find(lambda == cw_b),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L,Refmat_sol(find(lambda == cw_g),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L,Refmat_sol(find(lambda == cw_r),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity Under Solution')
xlabel('L (nm)','FontSize',16);
ylabel('Reflectivity ','FontSize',16');
saveas(figure(6),[pwd '/Figures/RefSim/6LvsRefSolRGB.fig']);

%% COMPARE THICKNESS VS REFLECTANCE UNDER WATER (n=1.33) AND SOLUTION (n=1.34)

figure(7)
hold on
plot(L,Refmat_water(find(lambda == cw_b),:),'b','LineWidth',2) %Water under B
plot(L,Refmat_sol(find(lambda == cw_b),:),'b--','LineWidth',2) %Solution under B
plot(L,Refmat_water(find(lambda == cw_g),:),'g','LineWidth',2) %Water under G
plot(L,Refmat_sol(find(lambda == cw_g),:),'g--','LineWidth',2) %Solution under G
plot(L,Refmat_water(find(lambda == cw_r),:),'r','LineWidth',2) %Water under R
plot(L,Refmat_sol(find(lambda == cw_r),:),'r--','LineWidth',2) %Solution under R
legend("Water (n = 1.33)","Solution (n = 1.34)")
title("Thickness vs Reflectance Curve for Water and Solution")
xlabel('L (nm)','FontSize',16);
ylabel('Reflectivity ','FontSize',16');
saveas(figure(7),[pwd '/Figures/RefSim/7LvsRefWater&SolRGB.fig']);

%% dR/dn AT 1.33-1.34 CHANGE IN REFLECTANCE AS n CHANGES FROM 1.33 TO 1.34 

dRdn_B = Refmat_water(find(lambda == cw_b),:) - Refmat_sol(find(lambda == cw_b),:); 
dRdn_G = Refmat_water(find(lambda == cw_g),:) - Refmat_sol(find(lambda == cw_g),:);
dRdn_R = Refmat_water(find(lambda == cw_r),:) - Refmat_sol(find(lambda == cw_r),:);

%% PLOT dR/dn (CHANGE IN REFLECTIVITY OVER AS n CHANGES FROM 1.33 TO 1.34)

figure(8)
hold on
plot(L,dRdn_B,'b','LineWidth',2)
plot(L,dRdn_G,'g','LineWidth',2)
plot(L,dRdn_R,'r','LineWidth',2)
xline(L(find(dRdn_R==min(dRdn_R))))
title("Thickness vs Change in Reflectance for n = 1.33-1.34")
xlabel('L (nm)','FontSize',16);
ylabel('dR / dn ','FontSize',16');
text(100,0.003,"Line = 108.6 nm")
saveas(figure(8),[pwd '/Figures/RefSim/8LvsdRdnRGB.fig']);

%% dR(t)/R (CHANGE IN REFLECTANCE AS THICKNESS CHANGES FROM 0 TO 300 nm)

diffL = diff(Refmat_air,[],2); %Reflectance change as L increases from 0 to 300 nm

dR_R_blue = diffL(lambda==cw_b,:); %Change in reflectance as L increases in blue
dR_R_green = diffL(lambda==cw_g,:); %Change in reflectance as L increases in green
dR_R_red = diffL(lambda==cw_r,:); %Change in reflectance as L increases in red

%% PLOT dR(t)/R
figure(9)
hold on
plot(L(2:end),dR_R_blue,'b','LineWidth',2) %blue
plot(L(2:end),dR_R_green,'g','LineWidth',2) %green
plot(L(2:end),dR_R_red,'r','LineWidth',2) %red
xlabel("d (nm)")
ylabel("Delta R(t) / R")
xlim([100 150])
title("Change in reflectance over change in thickness")
saveas(figure(9),[pwd '/Figures/RefSim/9LvsdRRRGB.fig']);

%% SPECTRUM DATA FOR RGB LEDS

load("Osram_Spec_Data_edited.mat")
bluespectrum = Spec_DentalBlue_460nm; clear("Spec_DentalBlue_460nm")
greenspectrum = Spec_Green_517nm; clear("Spec_Green_517nm")
redspectrum = Spec_Red_633nm; clear("Spec_Red_633nm")


%% INTERPOLATE SPECTRUM DATA

intbluespectrum = interp1(Osram_lambda,bluespectrum,lambda); 
intgreenspectrum = interp1(Osram_lambda,greenspectrum,lambda);
intredspectrum = interp1(Osram_lambda,redspectrum,lambda);

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
reflectivity_curve_0nm = Refmat_air(:,L==0); %Reflectivity curve at 0 nm (used below)
%%
for i = 1:numval
Ref_spec_red = Refmat_air(i,find(L==0))'.*intredspectrum; %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green = Refmat_air(i,find(L==0))'.*intgreenspectrum; %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue = Refmat_air(i,find(L==0))'.*intbluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end

%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)
figure(11)
hold on
plot(lambda,reflectivity_curve_0nm,'k-','LineWidth',2) %Reflectivity curve at 0 nm
plot(lambda,intredspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,intgreenspectrum,'g','Linewidth',2) %Green spectrum before multiplication
plot(lambda,intbluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication
plot(lambda,Ref_spec_red,'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green,'g','Linewidth',2) %Green spectrum after multiplication
plot(lambda,Ref_spec_blue,'b','Linewidth',2) %Blue spectrum after multiplication
xlim([400 700])
ylim([0 1.3])
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon (no oxide)')
legend('Reflectivity Curve for Silicon (0 nm)','Red Spectrum and I_O_u_t','Green Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','northeast')
saveas(figure(11),[pwd '/Figures/RefSim/11SiRef.fig']);

%% Calculate total reflected intensity by multiplying reflectance with reflected intensity (I_Out)
I_refred = []; %Reflected Intensity for red
I_refgreen = []; %Reflected Intensity for green
I_refblue = []; %Reflected Intensity for blue
for i=1:numval
   fprintf("Now running %.0f\n",i)
   I_refred(:,i) = sum(Refmat_air(2:end,i).*Ref_spec_red(2:end)');      
   I_refgreen(:,i) = sum(Refmat_air(2:end,i).*Ref_spec_green(2:end)');     
   I_refblue(:,i) = sum(Refmat_air(2:end,i).*Ref_spec_blue(2:end)');       
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

%% Display simulated color for 0-300 nm by comparing RGB reflected intensity values
figure(13)
I_tot = cat(3, I_refred, I_refgreen, I_refblue);
I_totnorm = (I_tot./max(I_tot)).*255;
for i = 1:6                               
   I_totnorm = [I_totnorm; I_totnorm]; %Add to each other for visualization purposes
end
imtot = uint8(round(I_totnorm));
imagesc(imtot)
set(gca,'YDir','normal')
ylim([1,30])
title('Simulated Color Si-SiO_2')
xlabel('Thickness')
xticks([find(L==0) find(L==50) find(L==100) find(L==150) find(L==200) find(L==250) find(L==300)])
xticklabels([0 50 100 150 200 25+0 300])
saveas(figure(13),[pwd '/Figures/RefSim/13SimColor.fig']);

%%

%% Variables needed for thuckness calc
Ref_Si= abs((n_Air(:,2) - n_Si(:,2).^2)./ (n_Air(:,2) + n_Si(:,2).^2)).^2 ;
%%
Refmat_B = Refmat_air/Ref_Si(find(lambda == cw_b));
Refmat_G = Refmat_air/Ref_Si(find(lambda == cw_g));
Refmat_R = Refmat_air/Ref_Si(find(lambda == cw_r));
%%
Ref_at_0nm_B = Refmat_air(find(lambda == cw_b),find(L == 0));
Ref_at_0nm_G = Refmat_air(find(lambda == cw_g),find(L == 0));
Ref_at_0nm_R = Refmat_air(find(lambda == cw_r),find(L == 0));

Ref_at_120nm_B = Refmat_air(find(lambda == cw_b),find(L == 120));
Ref_at_120nm_G = Refmat_air(find(lambda == cw_g),find(L == 120));
Ref_at_120nm_R = Refmat_air(find(lambda == cw_r),find(L == 120));

%%
save("SimOutputs/RefSimtoBulkData.mat",'lambda','L','cw_r','n_SiO2','n_Si','Refmat_air','I_refred','I_refgreen','I_refblue')
%%
save("SimOutputs/RefSimData")

%% Color for 120 nm
% figure(14)
% im_120 = imtot(:,find(L==120),:);
% for i = 1:6
%    im_120 = [im_120 im_120];
% end
% imagesc(im_120)

%TEST
%ALL2
%JAN11NIGHT




