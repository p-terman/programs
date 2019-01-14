function LIP_scratch(sim_list)
% this will get the energy per cm for events that have 5 s2s (of any type)
total = zeros(1000, length(sim_list));
ave = zeros(length(sim_list),1);
sd = zeros(length(sim_list),1);
ave_3SD = zeros(length(sim_list),1);
SD_3SD = zeros(length(sim_list),1);
ave_per_cm = zeros(length(sim_list),1);
sd_per_cm = zeros(length(sim_list),1);
figure 
hold on 

for i =1:length(sim_list)

    dp_name = sim_list{i};
    dp=load(dp_name);
    dp= which_class5_type (dp);
    simple_energy
    clear total cuts s2_cut total_cut
    total(:,i) = sum_total_E_keV;
    ave(i)=mean(total(:,i));
    sd(i)=std(total(:,i));
    total_holder = total(:,i);
    
  %  cut_size = total_holder<(ave(i)+sd(i)) & total_holder>(ave(i)-sd(i));
   % cut_nonZero = total_holder>0; - redundant 
    s2_pulse_cut = dp.x_corrected > -99;
    s2_cut = sum(s2_pulse_cut) >4; % adding in a five pulse requirement.
                                   %having a corrected pos means that it
                                   %met all previous requirements for expanded s2 s2
    cuts=s2_cut; %& cut_size';
    good_evts = find(cuts);
    total_cut = zeros (size(good_evts));
        for k=1:length(good_evts)
            j=good_evts(k);
            xx= dp.x_corrected(dp.x_corrected(:,j)>-99, j);
            yy= dp.y_corrected(dp.x_corrected(:,j)>-99, j);
            zz= dp.z_drift_samples(dp.x_corrected(:,j)>-99, j);
            S = find_length(xx, yy, zz);
            d(j,i) = S;
            total_cut(j) = total_holder(j)/d(j,i);
        end
    total_cut=total_cut(good_evts) ;   % this is per cm at this point
    %want to hist total_cut?
    total_cut = total_cut(~isnan(total_cut));
    
    ave_per_cm(i) = mean(total_cut);
    sd_per_cm(i) = std(total_cut);
    total_3SD_cut=total_cut(total_cut<(ave_per_cm(i)+3*sd_per_cm(i)) & total_cut>(ave_per_cm(i)-3*sd_per_cm(i)));
    ave_3SD(i) = mean(total_3SD_cut);
    SD_3SD(i)= std(total_3SD_cut);
    h = histogram(total_3SD_cut, 'BinWidth' , 1 , 'BinLimit',  [ max([ 0 (ave_3SD(i)-SD_3SD(i))]) (ave_3SD(i)+SD_3SD(i))] );
end

hold off

end
-
function Load_filtered_LIPs 
%This will load all the files of a type (charge) of sim
files = dir('*T0050.mat');
Energy_05.total=[];
TL_05.total=[];
eptl_05.total=[];
Chi2_05.total=[];

for i=1:length(files)
    i
    load(files(i).name);
    struct_name = ['s' files(i).name(11: end-4)];
    Energy_05.(struct_name)=outArray.sum_total_E_keV;
    Energy_05.total=[Energy_05.total ; Energy_05.(struct_name)(:)];
    TL_05.(struct_name) = outArray.track_lengths;
    TL_05.total=[TL_05.total ; TL_05.(struct_name)(:)];
    eptl_05.(struct_name) = outArray.sum_total_E_keV ./outArray.track_lengths;
    eptl_05.total=[eptl_05.total ; eptl_05.(struct_name)(:)];
    Chi2_05.(struct_name) = outArray.Chi2_dof_array;
    Chi2_05.total=[Chi2_05.total ; Chi2_05.(struct_name)(:)];

end
end

%load real data that has been processed with lip filter

files = dir('*.mat');
Energy_real.total=[];
TL_real.total=[];
eptl_real.total=[];
Chi2_real.total=[];
luxstamps.total = [];
nS2=[];

for i=1:length(files)
i
load(files(i).name);
struct_name = ['s' files(i).name(11: end-4)];
%  Energy_real.(struct_name)=outArray.sum_total_E_keV;
Energy_real.total=[Energy_real.total ; outArray.sum_total_E_keV'];
%  TL_real.(struct_name) = outArray.track_lengths;
TL_real.total=[TL_real.total ; outArray.track_lengths'];
% eptl_real.(struct_name) = outArray.sum_total_E_keV ./outArray.track_lengths;
eptl_real.total=[eptl_real.total ; (outArray.sum_total_E_keV ./ outArray.track_lengths)'];%  Chi2_real.(struct_name) = outArray.Chi2_dof_array;
Chi2_real.total=[Chi2_real.total ; outArray.Chi2_dof_array'];
luxstamps.total=[luxstamps.total ; outArray.luxstamps'];
nS2 = [nS2 ; outArray.nS2'];
end


%cuts to real

Chi2_real.cutChi2=[]
Chi2_real.cutChi2=Chi2_real.total(Chi2_real.total<1)
Chi2_real.cutTL=[]
Chi2_real.cutTL=Chi2_real.total(TL_real.total>40)
Chi2_real.cutTL=Chi2_real.total(TL_real.total>40 & Chi2_real.total<1)
Chi2_real.cutTLChi2=[]
Chi2_real.cutTLChi2=Chi2_real.total(TL_real.total>40 & Chi2_real.total<1)
Energy_real.cutTL=Energy_real.total(TL_real.total>40)
Energy_real.cutChi2=Energy_real.total(Chi2_real.total<1)
Energy_real.cutTLChi2=Energy_real.total(Chi2_real.total<1 & TL_real.total>40)
luxstamps_cut=[];
luxstamps_cut=luxstamps.total(TL_real.total>40 & Chi2_real.total<1);


%% looking at the cut real data and then putting it in a place by itself

tf = ismember(dp.luxstamp_samples, luxstamps_cut);
index = find(tf);
dp_cut.pulse_length_samples = dp.pulse_length_samples(:, index);
dp_cut.z_drift_samples = dp.z_drift_samples(:, index);
dp_cut.pulse_classification = dp.pulse_classification(:, index);
dp_cut.xyz_corrected_pulse_area_all_phe = dp.xyz_corrected_pulse_area_all_phe(:, index);
dp_cut.x_corrected = dp.x_corrected(:, index);
dp_cut.y_corrected = dp.y_corrected(:, index);
dp_cut.sd_phiXR = dp.sd_phiXR(:, index);
dp_cut.sd_radius_inf = dp.sd_radius_inf(:, index);
dp_cut.sd_radius_sup = dp.sd_radius_sup(:, index);
dp_cut.luxstamp_samples = dp.luxstamp_samples(:, index);
dp_cut.pulse_start_samples = dp.pulse_start_samples(:, index);
dp_cut.pulse_area_phe = dp.pulse_area_phe(:, index);
dp_cut.s1s2_pairing = dp.s1s2_pairing(:, index);
dp_cut.admin = dp.admin;
dp_cut.aft_t0_samples  = dp.aft_t0_samples(:, index);
dp_cut.event_number = dp.event_number(:, index);
dp_cut.top_bottom_ratio = dp.top_bottom_ratio(:, index);
dp_cut.x_cm=dp.x_cm(:, index);
dp_cut.y_cm = dp.y_cm(:, index);




%% graph selected events - not ready
x_int_cut =outArray.x_int_array(find(ismember(outArray.luxstamps, luxstamps_cut)));
x_slope_cut  = outArray.x_slope_array(find(ismember(outArray.luxstamps, luxstamps_cut)));
y_int_cut =outArray.y_int_array(find(ismember(outArray.luxstamps, luxstamps_cut)));
y_slope_cut  = outArray.y_slope_array(find(ismember(outArray.luxstamps, luxstamps_cut)));

figure
 hold on 
for i = 1:10
    ind = find(dp.luxstamp_samples== luxstamps_cut(i));
    cuts = dp.x_corrected(:,ind) > -99 & (dp.pulse_classification(:,ind) == 2 | dp.pulse_classification(:,ind) == 4 | dp.s2_like_class5(:,ind) ==1);
    eptl_cut = (dp.xyz_corrected_pulse_area_all_phe(:,ind) ./ dp.pulse_length_samples(:,ind)) > 1;
    all_cuts = cuts .* eptl_cut;
    all_cuts = all_cuts == 1;
    z=0:54;
    x=x_int_cut(i)+ x_slope_cut(i) * z; 
    y=y_int_cut(i)+ y_slope_cut(i) * z; 
    RGB = [randi(256) randi(256) randi(256)]/256 ;
    
   
    plot(x,z,'color' , RGB);
    scatter(dp.x_corrected(all_cuts, ind),  double(dp.z_drift_samples(all_cuts, ind)) * 1.51/10 *0.01 , 10 ,  RGB)
end
hold off



s1_cut = (dp.pulse_classification == 1 | dp.s1_like_class5) & dp.s1s2_pairing == 1;
s2_cut = (dp.pulse_classification == 2 | dp.pulse_classification == 4 | dp.s2_like_class5 == 1) & dp.x_corrected > -99;
s1_E_keV = dp.xyz_corrected_pulse_area_all_phe .* s1_cut * 13.7/1000 * 1/0.117;
s2_E_keV = dp.xyz_corrected_pulse_area_all_phe .* s2_cut * 13.7/1000 * 1/12.1;
en= s1_E_keV+s2_E_keV
tot_en = sum(en)
mean(tot_en)
%% ------------------------------------
data = struct_01.eptl.cutTLChi2;
ave = mean(data);
sd = std(data);
sd3_cut = data(data < (3*sd +ave) & data> (ave - 3*sd));
ave_3sd_cut = mean(sd3_cut)

%%
legend('MI' , '10 GeV')
title('Energy per Track Length with TL and Chi^2 Cuts')
xlabel('keV/cm')
ylabel('counts')


%% graph energies 


data = struct_25_h.eptl.cutTLChi2;
ave = mean(data);
sd = std(data);
sd3_cut = data(data < (3*sd +ave) & data> (ave - 3*sd));
ave_3sd_cut_25h = mean(sd3_cut)

data = struct_25.eptl.cutTLChi2;
ave = mean(data);
sd = std(data);
sd3_cut = data(data < (3*sd +ave) & data> (ave - 3*sd));
ave_3sd_cut20 = mean(sd3_cut)

%line([ave_3sd_cut_h, ave_3sd_cut_h], ylim, 'LineWidth', 2, 'Color', 'b');
line([ave_3sd_cut, ave_3sd_cut], ylim, 'LineWidth', 2, 'Color', 'r');

legend('25 GeV' , 'MI' , '3SD mean 25 GeV' , '3SD mean MI')
title('STD of S2 Pa vs Lfit/Ldata')
xlabel('Lfit/Ldata')
ylabel('std S2 PA')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

%log density plot?
set(gca,'colorscale','log')



%%%%%%%55555 draw cicle on a log graph
gscatter(log(all.first_last(ncuta)), log(all.stdS2PA(ncuta)), all.sim_or_live(ncuta) , 'rb ', '*.')
c = circles (0, 9, 2 , 'facecolor', 'none')
