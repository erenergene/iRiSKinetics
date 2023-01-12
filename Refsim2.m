clear;clc;
addpath(genpath(pwd))
%% LOAD RefSimData

load("RefSimData.mat")

%% Multiply LED spectrum with reflectivity curves to find reflected Intensity

reflectivity_curve_0nm = Refmat_air(:,L==0); %Reflectivity curve at 0 nm (used below)
%%
for i = 1:numval
Ref_spec_red = Refmat_air(i,L==0)'.*intredspectrum; %Reflectivity spectrum for red at 0 nm (R_r)
Ref_spec_green = Refmat_air(i,L==0)'.*intgreenspectrum; %Reflectivity spectrum for green at 0 nm (R_g)
Ref_spec_blue = Refmat_air(i,L==0)'.*intbluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end
%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)
figure(1)
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
title('Reflectivity Curve for Silicon at 0 nm (no oxide)')
legend('Reflectivity Curve for Silicon at 0 nm','Red Spectrum and I_O_u_t','Green Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','northeast')

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
figure(5)
hold on
plot(L,I_refblue/100,'b','LineWidth',2)
plot(L,I_refgreen/100,'g','LineWidth',2)
plot(L,I_refred/100,'r','LineWidth',2)
title('Reflected Intensity (I_O_u_t) for RGB from 0 nm to 300 nm')
legend('Blue','Green','Red')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);

%% Display simulated color for 0-300 nm by comparing RGB reflected intensity values
figure(6)
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
xticklabels([0 50 100 150 200 250 300])
%% Color for 120 nm
figure(7)
im_120 = imtot(:,find(L==120),:);
for i = 1:6
   im_120 = [im_120 im_120];
end
imagesc(im_120)