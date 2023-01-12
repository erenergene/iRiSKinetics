function est_L = reftocurve_lsqr(RefRed, RefGreen, RefBlue)

load("RefSimData.mat")

% Ref_at_120nm_B = 0.2868
% Ref_at_120nm_G = 0.1785
% Ref_at_120nm_R = 0.0866

%% Get reflectance points at RGB

RefBlue = [cw_b, RefBlue]; %Reflectance Point at Blue
RefGreen = [cw_g, RefGreen]; %Reflectance Point at Green
RefRed = [cw_r, RefRed]; %Reflectance Point at Red

%% Subtract reflectance value of theoretical value from the measured reflectance at RGB 
%% Creates a n-by-3 (3 = RGB) matrix with difference of measured and theoritical reflectance values. 

for i = 1:length(Refmat_air)
    Diff_at_B(:,i) = abs(RefBlue - Refmat_B(find(lambda==cw_b),i));
    Diff_at_G(:,i) = abs(RefGreen - Refmat_G(find(lambda==cw_g),i));
    Diff_at_R(:,i) = abs(RefRed - Refmat_R(find(lambda==cw_r),i));

%     Diff_at_B(:,i) = abs(RefBlue - Gamma(valtoindex(cw_b,numval,Osram_lambda(1),Osram_lambda(end)),i));
%     Diff_at_G(:,i) = abs(RefGreen - Gamma(valtoindex(cw_g,numval,Osram_lambda(1),Osram_lambda(end)),i));
%     Diff_at_R(:,i) = abs(RefRed - Gamma(valtoindex(cw_r,numval,Osram_lambda(1),Osram_lambda(end)),i));

    L_Diff(i,:) = [Diff_at_B(2,i) Diff_at_G(2,i) Diff_at_R(2,i)]; 

end

lsqr_sum_L_Diff = sum(((L_Diff).^2),2); %average of R,G,B differences at each thickness
est_in = find(lsqr_sum_L_Diff == min(lsqr_sum_L_Diff)); %find index of minimum difference
est_L = L(est_in); %find estimated thickness

end