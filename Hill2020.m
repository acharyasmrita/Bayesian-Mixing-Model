% Supplementary data use from Hill et al., 2020


fid = fopen ('Hill20.csv');
readData = textscan (fid,'%f %f %f %f %f %f %f %f %f %f', 'HeaderLines',1,'Delimiter',',');
vals = csvread ('Hill20.csv',1,0); day = vals (:,1);
discharge = vals (:,2);
DO = vals(:,3);
DH = vals(:,4);
O_ice = vals (:,5);
H_ice = vals (:,6);
O_rain = vals (:,7);
H_rain = vals (:,8);
O_snow = vals (:,9);
H_snow = vals (:,10);
sigmaO = 1.8;
sigmaH = 15.3;
N = 2e7;
sample_number = 1:size(vals,1);
sigma_fiO = 0.2;
sigma_fiH = 1.6;
sigma_frO = 4.1;
sigma_frH = 33.5;
sigma_fsO = 5.2;
sigma_fsH = 39.8;
fi_best = zeros(size(discharge));
fr_best = zeros(size(discharge));
fs_best = zeros(size(discharge));
fi_mean = zeros(size(discharge));
fr_mean = zeros(size(discharge));
fs_mean = zeros(size(discharge));
fi_std = zeros(size(discharge));
fr_std = zeros(size(discharge));
fs_std = zeros(size(discharge));
fifr_corr = zeros(size(discharge));
frfs_corr = zeros(size(discharge));
fsfi_corr = zeros(size(discharge));
DOp_best = zeros(size(discharge));
DHp_best = zeros(size(discharge));
for n = sample_number
    fs_trial = unifrnd(0,1,N,1);
    fi_trial = unifrnd(0,1,N,1);
    pos = find((fs_trial + fi_trial)<=1);
    fs_trial = fs_trial(pos);
    fi_trial = fi_trial(pos);
    fr_trial = 1-(fs_trial + fi_trial);
    O_ice_trial = normrnd(mean(O_ice),sigma_fiO,length(pos),1);
    O_rain_trial = normrnd(mean(O_rain),sigma_frO,length(pos),1);
    O_snow_trial = normrnd(mean(O_snow),sigma_fsO,length(pos),1);
    DOp = (fs_trial.* O_snow_trial)+(fi_trial.* O_ice_trial)+(fr_trial.* O_rain_trial);
    H_ice_trial = normrnd(mean(H_ice),sigma_fiH,length(pos),1);
    H_rain_trial = normrnd(mean(H_rain),sigma_frH,length(pos),1);
    H_snow_trial = normrnd(mean(H_snow),sigma_fsH,length(pos),1);
    DHp = (fs_trial.* H_snow_trial)+(fi_trial.* H_ice_trial)+(fr_trial.* H_rain_trial);
    phi = (DOp - DO (n)).^2./sigmaO.^2 + (DHp - DH (n)).^2./sigmaH.^2;
    [best_phi,pos] = min(phi);
    fi_best(n) = fi_trial(pos);
    fr_best(n) = fr_trial(pos);
    fs_best(n) = fs_trial(pos);
    DOp_best(n)= DOp(pos);
    DHp_best(n)= DHp(pos);
    L = exp((-1/2).*phi); L=L./max(L);
    accept=find(L>rand(size(L)));
    fprintf(1,['day %d:tried %d samples of the prior,accepted %d','samples of the posterior; %6.2f%s acceptance rate/n'],n,length(fs_trial),length(accept),100*length(accept)/length(fs_trial),char(37));
    phi_post = phi(accept);
    fi_post = fi_trial(accept);
    fr_post = fr_trial(accept);
    fs_post = fs_trial(accept);
    O_ice_post = O_ice_trial(accept);
    H_ice_post = H_ice_trial (accept);
    O_rain_post = O_rain_trial(accept);
    H_rain_post = H_rain_trial(accept);
    O_snow_post = O_snow_trial(accept);
    H_snow_post = H_snow_trial(accept);
    save(sprintf('Post_%03d.mat',n),'phi_post','fi_post','fr_post','fs_post','O_ice_post','H_ice_post','O_rain_post','H_rain_post','O_snow_post','H_snow_post');
    fi_mean(n)= mean(fi_trial(accept));
    fr_mean(n)= mean(fr_trial(accept));
    fs_mean(n)= mean(fs_trial(accept));
 Cmat = cov([fi_trial(accept),fr_trial(accept),fs_trial(accept)]);
 fi_std(n)=sqrt(Cmat(1,1));
 fr_std(n)=sqrt(Cmat(2,2));
 fs_std(n)=sqrt(Cmat(3,3));
 fifr_corr(n)= Cmat(1,2)/(fi_std(n)*fr_std(n));
 fifs_corr(n)= Cmat(1,3)/(fi_std(n)*fs_std(n));
 frfs_corr(n)= Cmat(2,3)/(fr_std(n)*fs_std(n));
end
clear fi_trial fr_trial fs_trial;
clear O_ice_trial H_ice_trial;
clear O_rain_trial H_rain_trial;
clear L accept phi;
clear phi_post fs_post fi_post fr_post;
clear O_snow_post H_snow_post;
clear O_rain_post H_rian_post;
discharge_i = fi_best.*discharge;
discharge_r = fr_best.*discharge;
discharge_s = fs_best.*discharge;
save OHL_soln.mat
% Return

post_1 = load('Post_001.mat');
post_2 = load('Post_002.mat');
solution_mat = load('OHL_soln.mat');

figure; 
h1 = plot_gaussian_ellipsoid([ solution_mat.fi_mean(1,1) solution_mat.fr_mean(1,1)], [Cmat(2, 2:3); Cmat(3, 2:3)]);
h2 = plot_gaussian_ellipsoid([ solution_mat.fi_mean(2,1) solution_mat.fr_mean(2,1)], [Cmat(2, 2:3); Cmat(3, 2:3)]);
set(h1,'color','r');
set(h2, 'color','b');
xlabel('fr');
ylabel('fi');
title('Day 1 and Day 2 Fi vs Fr Gaussian Ellipses');
 
figure; 
h1 = plot_gaussian_ellipsoid([ solution_mat.fr_mean(1,1) solution_mat.fs_mean(1,1)], [Cmat(2, 2:3); Cmat(3, 2:3)]);
h2 = plot_gaussian_ellipsoid([ solution_mat.fr_mean(2,1) solution_mat.fs_mean(2,1)], [Cmat(2, 2:3); Cmat(3, 2:3)]);
set(h1,'color','r');
set(h2, 'color', 'b');
xlabel('fr');
ylabel('fs');
title('Day 1 and Day 2 Fi vs Fr Gaussian Ellipses');



