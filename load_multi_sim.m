function dp = load_multi_sim(sim)

     
%% load rqs 
 
dataset = {sim};

%%
for jj = 1:length(dataset)
   % datapath=[dataset{jj} filesep 'matfiles']
        %datapath = strcat('/eliza3/lux/rq/',dataset{jj} );
      % datapath = ['/media/paul/TOSHIBA/100_pulse/' filesep dataset{jj} filesep 'matfiles'] 
        rqm_path = strsplit(sim, '/');; 
    %     rqs = {'all'};
 
         rqs = {};
         rqs{end+1} = 'pulse_area_phe';
        %  rqs{end+1} = 'peak_area_phe';
        % rqs{end+1} = 'file_number';
         rqs{end+1} = 'event_number';
         rqs{end+1} = 'admin';
         rqs{end+1} = 'pulse_length_samples';
       %  rqs{end+1} = 'n_samples_in_evt';
        rqs{end+1} = 'z_drift_samples';
         rqs{end+1} = 'pulse_classification';
         rqs{end+1} = 'xyz_corrected_pulse_area_all_phe';
         rqs{end+1} = 'xyz_corrected_pulse_area_bot_phe';
         rqs{end+1} = 'x_corrected';
         rqs{end+1} = 'y_corrected';
         rqs{end+1} = 'x_cm';
         rqs{end+1} = 'y_cm';
         rqs{end+1} = 's1s2_pairing';
        %rqs{end+1} = 'chi2';
     %    rqs{end+1} = 'rec_dof';
         rqs{end+1} = 'pulse_length_samples';
        %rqs{end+1} = 'event_timestamp_samples';
         %rqs{end+1} = 'peak_height_mV';
         %rqs{end+1} = 'hft_t10l_samples';
         %rqs{end+1} = 'hft_t10r_samples';
       %  rqs{end+1} = 's2filter_max_area_diff';
        rqs{end+1} = 'aft_t0_samples';
       % rqs{end+1} = 'aft_t1_samples';
       % rqs{end+1} = 'aft_t2_samples';
     %   rqs{end+1} = 'aft_t25_samples';
      %  rqs{end+1} = 'aft_t75_samples';
     %  rqs{end+1} = 'exp_fit_tau_fall_samples';
        %rqs{end+1} = 'livetime_latch_samples';
       % rqs{end+1} = 'livetime_end_samples';
      % rqs{end+1} = 'full_evt_area_phe';
       rqs{end+1} = 'sd_phiXR';
       rqs{end+1} = 'sd_radius_inf'; 
       rqs{end+1} = 'sd_radius_sup'; 
        rqs{end+1} = 'luxstamp_samples';
         rqs{end+1} = 'pulse_start_samples';
       % rqs{end+1} = 'top_bottom_ratio';
       % rqs{end+1}= 'prompt_fraction_tlx';
        rqs{end+1} = 'top_bottom_ratio';
       % rqs{end+1} = 'peak_height_phe_per_sample';
       % rqs{end+1} = 'skinny_peak_area_phe';
       % rqs{end+1} = 'prompt_fraction';
       % rqs{end+1} = 'pulse_height_phe_per_sample';
       % rqs{end+1} = 'rms_width_samples';
       % rqs{end+1} = 'selected_s1_s2';
       % rqs{end+1} = 'multiple';
        %rqs{end+1} = 'spike_count';
        %rqs{end+1} = 'correction_electron_lifetime';
        %rqs{end+1} = 'correction_s2_xy_dependence';
        %rqs{end+1} = 'correction_s1_xyz_dependence';
        %rqs{end+1} = 'correction_s2_xy_dependence_bot';
        % rqs{end+1} = 'correction_s1_xyz_dependence_bot';
        %rqs{end+1} = 'event_with_dd_trigger';
        %rqs{end+1} = 'dd_trigger_pulse_area';
        %rqs{end+1} = 'dd_pulse_start_sample';
        % luxstamp_samplesrqs{end+1} = 's1_like_class5';
       %   rqs{end+1} = 's1_like_class5';
        %  rqs{end+1} = 's2_like_class5';
       
        num_files_max = [2000];
        dp_temp = LUXLoadMultipleRQMs_framework(rqm_path,rqs,num_files_max);    
        
        
        % apply position reconstruction corrections
        
%         dp_temp = Corrections_PositionCorrection_Function(dp_temp)
        
        % change the rq variables type
        
%         dp_temp.dataset_number = uint16( use_case .* ones(1,length(dp_temp.event_number)) );
        dp_temp.pulse_area_phe = single(dp_temp.pulse_area_phe);
        dp_temp.event_number = single(dp_temp.event_number);
        dp_temp.pulse_classification = uint8(dp_temp.pulse_classification);
        dp_temp.x_corrected = single(dp_temp.x_corrected);
        dp_temp.y_corrected = single(dp_temp.y_corrected);
        dp_temp.x_cm = single(dp_temp.x_cm);
        dp_temp.y_cm = single(dp_temp.y_cm);
        dp_temp.xyz_corrected_pulse_area_all_phe = single(dp_temp.xyz_corrected_pulse_area_all_phe);
        dp_temp.xyz_corrected_pulse_area_bot_phe = single(dp_temp.xyz_corrected_pulse_area_bot_phe);
       % dp_temp.event_timestamp_samples = single(dp_temp.event_timestamp_samples);
        dp_temp.aft_t0_samples = single(dp_temp.aft_t0_samples);
        %dp_temp.aft_t2_samples = single(dp_temp.aft_t2_samples);
        %dp_temp.top_bottom_ratio = single(dp_temp.top_bottom_ratio);
        %dp_temp.dataset_number = uint16(dp_temp.dataset_number);
       % dp_temp.exp_fit_tau_fall_samples = single(dp_temp.exp_fit_tau_fall_samples);
       % dp_temp.golden = uint8(dp_temp.golden);
       
       % dp_temp.rms_width_samples = single(dp_temp.rms_width_samples);
        dp_temp.z_drift_samples = single(dp_temp.z_drift_samples);
        dp_temp.s1s2_pairing = single(dp_temp.s1s2_pairing);
       % dp_temp.selected_s1_s2 = uint8(dp_temp.selected_s1_s2);
        %dp_temp.multiple = uint8(dp_temp.multiple);  
        %dp_temp.spike_count = single(dp_temp.spike_count)
        %dp_temp.correction_electron_lifetime = single(dp_temp.correction_electron_lifetime);
        %dp_temp.correction_s2_xy_dependence = single(dp_temp.correction_s2_xy_dependence);    
        %dp_temp.correction_s1_xyz_dependence = single(dp_temp.correction_s1_xyz_dependence);
        %dp_temp.correction_s2_xy_dependence_bot = single(dp_temp.correction_s2_xy_dependence_bot);    
        % dp_temp.correction_s1_xyz_dependence_bot = single(dp_temp.correction_s1_xyz_dependence_bot);
           
        
        % clear the large peak_area_phe and spike_count rq or change types
    %{
        load_per_ch = true;
        if ~load_per_ch
                % dp_temp = rmfield(dp_temp,'peak_area_phe');
                dp_temp = rmfield(dp_temp,'spike_count');
        end
    
        if load_per_ch
            dp_temp.spike_count = uint32(dp_temp.spike_count);
            % dp_temp.peak_area_phe = single(dp_temp.peak_area_phe);
        end
      %}  
        if exist('dp','var')
           dp = LUXConcatdp(dp,dp_temp)
           dp.livetime_latch_samples = double(dp.livetime_latch_samples);
           dp.livetime_end_samples = double(dp.livetime_end_samples);
           lt = sum(dp.livetime_end_samples-dp.livetime_latch_samples)*10/1e9; %livetime in seconds
           lt_hour = lt/3600
        else
           dp = dp_temp
           
        end
        clear dp_temp;
 
        % save('Run4DDdp11205','dp','-v7.3')
end
 
 
% save(sprintf('Run4DDmultiz7p5dp%d',ii),'dp','-v7.3')
% dp = Corrections_PositionCorrection_Function(dp)
% save('20131218TritiumDP2p0VUVspike','dp','-v7.3')
 
dp.livetime_latch_samples = double(dp.livetime_latch_samples);
dp.livetime_end_samples = double(dp.livetime_end_samples);
lt = sum(dp.livetime_end_samples-dp.livetime_latch_samples)*10/1e9; %livetime in seconds
lt_hour = lt/3600;
 
if jj == 1
    total_dp = dp;
else
    total_dp=LUXConcatdp(total_dp,dp);
end

dp=total_dp;
%{
%%S1 energy calculation
energy_per_quanta=13.7; %eV
efficency=.12; %As a percent, amount fo light collected
s1_energy = dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S1  & cut_e_drift_time & cut_e_radius)  ./ efficency .* energy_per_quanta;
 
%%S2 energy calculation
Extraction_eff= 0.7 ;% eff of light collected 
phe_per_electron=27;
s2_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S2 & cut_e_drift_time & cut_e_radius) ./Extraction_eff .*energy_per_quanta ./phe_per_electron;
 
%%sphe energy
sphe_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_sphe & cut_e_drift_time & cut_e_radius) .*energy_per_quanta ./efficency;
 
%%single_e energy
single_e_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_single_e & cut_e_drift_time & cut_e_radius).*energy_per_quanta ./phe_per_electron;
 
total_energy = s1_energy + s2_energy + sphe_energy + single_e_energy;

%%S1 energy calculation
energy_per_quanta=13.7; %eV
efficency=.12; %As a percent, amount fo light collected
s1_energy = dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S1  & cut_e_drift_time & cut_e_radius)  ./ efficency .* energy_per_quanta;
 
%%S2 energy calculation
Extraction_eff= 0.7 ;% eff of light collected 
phe_per_electron=27;
s2_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S2 & cut_e_drift_time & cut_e_radius) ./Extraction_eff .*energy_per_quanta ./phe_per_electron;
 
%%sphe energy
sphe_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_sphe & cut_e_drift_time & cut_e_radius) .*energy_per_quanta ./efficency;
 
%%single_e energy
single_e_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_single_e & cut_e_drift_time & cut_e_radius).*energy_per_quanta ./phe_per_electron;
 
total_energy = s1_energy + s2_energy + sphe_energy + single_e_energy;
 
%%S1 energy calculation
energy_per_quanta=13.7; %eV
efficency=.12; %As a percent, amount fo light collected
s1_energy = dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S1  & cut_e_drift_time & cut_e_radius)  ./ efficency .* energy_per_quanta;
 
%%S2 energy calculation
Extraction_eff= 0.7 ;% eff of light collected 
phe_per_electron=27;
s2_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_golden_single_scatter_S2 & cut_e_drift_time & cut_e_radius) ./Extraction_eff .*energy_per_quanta ./phe_per_electron;
 
%%sphe energy
sphe_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_sphe & cut_e_drift_time & cut_e_radius) .*energy_per_quanta ./efficency;
 
%%single_e energy
single_e_energy= dp.xyz_corrected_pulse_area_all_phe .* (cut_p_single_e & cut_e_drift_time & cut_e_radius).*energy_per_quanta ./phe_per_electron;
 
total_energy = s1_energy + s2_energy + sphe_energy + single_e_energy;
 
summed_energy= sum(total_energy);
%%plot it
histogram(summed_energy)
set(gca, 'yscale', 'log')
%}