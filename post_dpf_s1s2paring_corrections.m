filename_evt='luxsm_20130624T0010_f000000001.evt'%'luxsm_20171016T1938'
%this name doesn't really need to bve exact, just to get right cal data
filename_prefix = filename_evt(1:19);
addpath(genpath(['~/LUXcode/Stable_Releases/DataProcessingStableReleases/v2.0']));
pathbase= '/home/paul/LUXdata' 
data_path_evt = [pathbase filesep filename_prefix filesep];
filename_rq = strrep(filename_evt,'evt','rq');
data_path_rq = [pathbase filesep];


fprintf('\n\n *** check dp_settings_new for max num pulses ******* \n');

 lug_iqs_xml_file = which('lug_iqs_new4.xml');
iq_xml_path=lug_iqs_xml_file;

data_processing_settings_path = which('dp_settings_new.xml');
data_processing_xml_path =data_processing_settings_path;

%dp1 = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

if isfield(dp, 'livetime_latch_samples');
    livetime=struct('livetime_latch_samples' , dp.livetime_latch_samples , 'livetime_end_samples' ,  dp.livetime_end_samples); 
    dp.admin.livetime= livetime;
else
    livetime=struct('livetime_latch_samples' , dp.admin.livetime.livetime_latch_samples , 'livetime_end_samples' ,  dp.admin.livetime.livetime_end_samples); 
end

%% Read the LUG settings and the table with the corrections


dp_settings_xml = XMLReader_framework(data_processing_xml_path); %load dp_settings
lug_iqs_xml = XMLReader_framework(iq_xml_path); %load lug iqs
settings.evt_settings = dp.admin.evt_settings; %grabs the settings arrays and makes a new array out of that
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;


%% S1S2 pairing

myname = 'S1S2Pairing_Naive';
module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

% Initialize variables


N = length(dp.event_number);
pulse_event_size = [max_num_pulses N];

dp.s1s2_pairing    = zeros(pulse_event_size, 'uint32');
dp.z_drift_samples =   zeros(pulse_event_size, 'uint32');


dp = which_class5_type(dp);  %180205 PAT added 

%% Loop per event and assign S1 S2 pairs

for ii = 1:N
    
    % find the indices for S1 and S2 pulses
    s1_inds = find(dp.pulse_classification(:,ii) == 1); %180205 only s1 here | dp.pulse_classification(:,ii) == 5); %treating all class 5 as an s1 -  PAT
    if isempty(s1_inds) % no normal s1
        s1_inds = find(dp.s1_like_class5(:, ii)==1); %use class5 as s1 if there
    end
    s2_inds = find(dp.pulse_classification(:,ii) == 2 | dp.pulse_classification(:,ii) == 4 | dp.s2_like_class5(:,ii)==1); % treating all SE as s2. 180205 added s2_like_class5 as well

    
    % how many pulses of each type are there?
    num_s1 = length(s1_inds);
    num_s2 = length(s2_inds);
   
    % are there one or more pulses of each type?
    has_s1 = num_s1 > 0; 
    has_s2 = num_s2 > 0;
   
    
    
    if has_s1 && has_s2  % if an s1/other and s2/SE is extant, pair the s1 to the s2s
        
        % find the first s1
        first_s1_ind = find(dp.pulse_classification(:,ii) == 1 | dp.s1_like_class5(:,ii)==1 ,1); %180205 brought back and added s1 like class5
             s2_after_s1 = find(dp.pulse_classification(:,ii) == 2 & (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ...   % after the first s1  - this is using aft_t0 as the 'start' time of the pulse - thus this just means that the s1 is timed before the s2           
            | dp.pulse_classification(:,ii) == 4 & (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ...  % now allows for SE
            | dp.s2_like_class5(:,ii)==1 &   (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ); %180205 now includes s2 like class5 pulses (logical structure may be redundent but I preserved for consistancy)       
        %takes all s2s following 'THE s1' 
        num_s2_after_s1 = length(s2_after_s1);

        if num_s2_after_s1 > 0 %are there s2s that came AFTER the s1

            first_s2_ind = s2_after_s1(1); %this is the first s2 immediately following the s1 (even if other types of pulses inbetween)


            % we choose a single s1 to pair with all s2 pulses (here, the first s1)
            dp.s1s2_pairing(first_s1_ind,ii) = 1;
            
            % loop over all s2 after the first s1, pairing them with the s1, and calculating drift time
            for ss = 1:num_s2_after_s1
                
                % label as a member of the pair - that this s2 has a
                % corresponding s1
                dp.s1s2_pairing(s2_after_s1(ss),ii) = 1;
                
                % calculate a drift time with respect to the s1
                delta_t = dp.aft_t0_samples(s2_after_s1(ss),ii) - dp.aft_t0_samples(first_s1_ind,ii);
                dp.z_drift_samples(s2_after_s1(ss),ii) = delta_t; %save drift time info
                  
            end
            
            
        end
        
    end
    
end


%% Corrections_PositionCorrection

myname = 'Corrections_PositionCorrection';
for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'xy_rec_cor')
            lrf_iq = lug_iqs_xml.iq(ii);
        end
    end
end %finds iq_type for the xy_rec_cor correction struct, tells its location and copies associated struct to lfr_iq

if ~isstruct(lrf_iq)
    disp(fprintf('\n\n %s: The IQ xy_rec_cor was not found\n\n************* Fatal Error*************\n\n',myname));
    return
end

position_correction_path = which('Corrections_PositionCorrection'); % find file location of Corrections_PositionCorrection
 
if ~isfield(lrf_iq, 'file_table_Mercury') & isfield(lrf_iq, 'file_table') % naming check, might not need
    lrf_iq.file_table_Mercury = lrf_iq.file_table;
end

if isfield(lrf_iq, 'file_table_Mercury') %this recovers the location of the cor_xy_iq correction array (it is a .mat file)
    position_correction_table_path = [position_correction_path(1:(end-numel(myname)-2)) lrf_iq.file_table_Mercury];
    if ~any(position_correction_table_path)
        disp(fprintf('\n\n %s: The file %s was not found in %s\n\n*************\n\n',myname, lrf_iq.file_table, position_correction_path(1:(end-numel(myname)-2))));
    else
        table_corrections = load(position_correction_table_path); % load the .mat corrections file
        %% Run the function. 
        dp = Corrections_PositionCorrection_Function(dp, table_corrections); %runs the main function
        
    end
    dp.x_corrected = single(dp.x_corrected); %type conversion
    dp.y_corrected = single(dp.y_corrected);

end



%% CorrectionsApplyCorrections

myname = 'ApplyCorrections';

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;
pmt_chs = 1:122;


[a1 b1] = size(dp.pulse_area_phe);  
[a2 b2] = size(dp.top_bottom_ratio);  adif = a1-a2;

blanks = [];

for i=1:b1
   blanks = [blanks 0]; 
end

if adif>0
   for i=1:adif
     dp.top_bottom_ratio = vertcat(dp.top_bottom_ratio,blanks);    
   end
end

% Initializing corrected values of area

dp.z_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.xyz_corrected_pulse_area_all_phe = dp.pulse_area_phe;

dp.z_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);
dp.xyz_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);

dp.correction_electron_lifetime = zeros(max_num_pulses,N, 'single');
dp.correction_s1_z_dependence = zeros(max_num_pulses,N, 'single');
dp.correction_s1_xyz_dependence = zeros(max_num_pulses,N, 'single');
dp.correction_s2_xy_dependence = zeros(max_num_pulses,N, 'single');

s1_ref_z_ns = mymodule_settings.detector_centre_z_ns;
allowed_gs = mymodule_settings.allowed_gs;


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Reading iq values from the IQs xml ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

[a num_iqs] = size(lug_iqs_xml.iq);

%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

lifetime_values = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'electron_attenuation_us').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
          lifetime_values = [lifetime_values lug_iqs_xml.iq(i).correction.fit.electron_attenuation_us];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
              
    end
  end  
 end
 
current_data_set_time=filename2epoch_framework(filename_prefix );
 
 if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
     [electron_lifetime] = InterpIQ(filename_prefix,dataset_times,lifetime_values);
 else    
     [index electron_lifetime_time] = NearestIQ(filename_prefix,dataset_times);
     electron_lifetime = lifetime_values(index);
 end
 
%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = zeros(0,3);
z_dep_bottom_values = zeros(0,3);
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'s1_both_zdep_quad_fit').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
           z_dep_both_values = vertcat(z_dep_both_values,lug_iqs_xml.iq(i).correction.fit.s1_both_zdep_quad_fit);
           z_dep_bottom_values = vertcat(z_dep_bottom_values,lug_iqs_xml.iq(i).correction.fit.s1_bottom_zdep_quad_fit);
           
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
            
    end
  end  
 end
 
   [index iq_time_zDep] = NearestIQ(filename_prefix,dataset_times);
   
   z_dep_par_all = z_dep_both_values(index,:);
   z_dep_par_bot = z_dep_bottom_values(index,:);
      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

s2_xy_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s2_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
           s2_xy_index = [s2_xy_index i];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];              
    end
  end  
 end
 
   [index iq_time_s2xyDep] = NearestIQ(filename_prefix,dataset_times);
      
   s2_x_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.x_bin_center;
   s2_y_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.y_bin_center;
   s2_map_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_both;
   s2_map_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_bottom;
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------
   
s1_xy_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
           s1_xy_index = [s1_xy_index i];                      
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];          
    end
  end  
 end
 
   [index iq_time] = NearestIQ(filename_prefix,dataset_times);
      
   s1_x_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.x_bin_center;
   s1_y_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.y_bin_center;
   s1_map_all = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_both;
   s1_map_bottom = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_bottom;

%------------------------------------------
% Finding the S1 xyz correction map values
%------------------------------------------
   
s1_xyz_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both_xyz').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
 %      if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
 %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);

           s1_xyz_index = [s1_xyz_index i];                      
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
%       else
%           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
%       end
              
    end
  end  
 end
 
   [index iq_time_s1xyzDep] = NearestIQ(filename_prefix,dataset_times);
      
   s1_xyz_x_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.x_bin_center;
   s1_xyz_y_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.y_bin_center;
   s1_xyz_z_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.z_bin_center;
   
   xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
   
   s1_xyz_map_all = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
   s1_xyz_map_bottom = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
    


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Calculating corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

%% Finding the drift time, x and y associated with S1

% Finding the drift time and xy position from the largest S2 pulse in the
% pairing


[a b] = size(dp.z_drift_samples); % b is number of events in file ??

drift = dp.z_drift_samples;
drift(find(isnan(drift))) = 0.0;
s1_drift_ns = +(dp.pulse_classification==1 | dp.s1_like_class5==1); % changed from only class 1 PAT 180131 - 180205 changed to have s1 and d2 like class 5
s1_x_cm = +(dp.pulse_classification==1 | dp.s1_like_class5==1); % same
s1_y_cm = +(dp.pulse_classification==1 | dp.s1_like_class5==1); %same
s2_phe =  dp.pulse_area_phe.*(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5==1); % changed from only class 2
s2_phe(find(isnan(s2_phe))) = 0.0;


 s1s = (sum(dp.s1s2_pairing==1)>1); % tells you if its an event with at least one s1 and one s2.
  
 for i=1:length(s1s)  % for each event in the file
     if s1s(i)>0 % if the event had an S1
       if s1s(i)==1   
        
         [v r] = max(s2_phe(:,i));% value and index (r for row) of Max S2
         [c1 r1 v1] = find(dp.pulse_classification(:,i)==1 | dp.s1_like_class5(:,1),1,'first');
   
         s1_drift_ns(c1,i) = 10.*drift(r,i); %row c1,r and column i
         s1_x_cm(c1,i) = dp.x_cm(r,i); %it seems that this is taking the location of the max s2 pulse and taking that to be the location of the s1
         %this can make sense since caluclated s1 pos is useless, so we
         %need a basis in the s2 to have a pos for s1
         s1_y_cm(c1,i) = dp.y_cm(r,i);  
         
       else %no s1 or class 5
          s1_drift_ns(:,i) = 0;
          s1_x_cm(:,i) = 0;
          s1_y_cm(:,i) = 0;
       end  
     end
 end

 
 
fprintf('Done finding S1 depths :) \n');

s1_drift_ns(find((s1_drift_ns==0))) = nan;
s1_x_cm(find((s1_x_cm==0))) = nan;
s1_y_cm(find((s1_y_cm==0))) = nan;



%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction

s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;

% Calculate electron lifetime correction (S2 Z-correction)
% Reading the values of electron lifetime from the LUG and interpolating
% between the two either side

electron_lifetime_correction = exp(double(dp.z_drift_samples)./(100.*electron_lifetime));
electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
 
%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

% Calculate S1 xyz corrections from map stored in IQ

s1xyz_correction = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_all,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction(find(isnan(s1xyz_correction))) = 1.0;

s1xyz_correction_bot = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_bottom,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction_bot(find(isnan(s1xyz_correction_bot))) = 1.0;

% Calculate S2 XY corrections

s2x = dp.x_cm.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5==1)); s2x(find(s2x==0)) = nan; %PAT changed to allow for class 4 to be counted as s2, and also s2 like class 5
s2y = dp.y_cm.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5==1)); s2y(find(s2y==0)) = nan; %180131 same - these seem to not be used??

s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'spline');
s2xy_correction = s2xy_correction.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5==1)); % PAT changed
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'spline');
s2xy_correction_bot = s2xy_correction_bot.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5==1));% PAT modified
s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;



% add RQs of correction factors

dp.correction_electron_lifetime = single(electron_lifetime_correction);
dp.correction_s1_z_dependence = single(s1_z_correction);
dp.correction_s1_z_dependence_bot = single(s1_z_correction_bot);
dp.correction_s1_xyz_dependence = single(s1xyz_correction);
dp.correction_s1_xyz_dependence_bot = single(s1xyz_correction_bot);
dp.correction_s2_xy_dependence = single(s2xy_correction);
dp.correction_s2_xy_dependence_bot = single(s2xy_correction_bot);

dp.admin.corrections.electron_lifetime = electron_lifetime;
dp.admin.corrections.s1_z_dependence.iq_time = iq_time_zDep;
dp.admin.corrections.s1_z_dependence.all = z_dep_par_all;
dp.admin.corrections.s1_z_dependence.bottom = z_dep_par_bot;
dp.admin.corrections.s1_xyz_dependence.iq_time = iq_time_s1xyzDep;
dp.admin.corrections.s2_xy_dependence.iq_time = iq_time_s2xyDep;


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Applying corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

dp.z_corrected_pulse_area_all_phe = dp.z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s1_z_correction; 
dp.z_corrected_pulse_area_bot_phe = dp.z_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s1_z_correction_bot; 

%--------------------------------------------------------------------------
% Apply XYZ corrections
%--------------------------------------------------------------------------

dp.xyz_corrected_pulse_area_all_phe = dp.xyz_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xyz_correction;
dp.xyz_corrected_pulse_area_bot_phe = dp.xyz_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s2xy_correction_bot.*s1xyz_correction_bot;

dp.s1_like_class5 = double(dp.s1_like_class5);
dp.s2_like_class5 = double(dp.s2_like_class5);

fprintf('Done!\n')

%status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime);

