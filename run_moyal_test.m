
input_path_filtered = '/home/paul/matfiles/dp_plus_chop/filtered/'
files = dir([input_path_filtered filesep 'lux*.mat']);
input_path_rq = '/home/paul/matfiles/dp_plus_chop/'

outputpath = [input_path_filtered filesep 'after_moyal/'];
for j = 1:length(files)

    current_file = files(j).name
    load([input_path_filtered filesep current_file]);
    current_file = current_file(1:end-4);
%    a = strcat(current_file, '_f000000001.rq.mat');
    %data = load(['/home/paul/LUXdata/' a]);
    %  list = dir([input_path_rq current_file(1:end-4) '*'])
    find_it = dir([input_path_rq current_file '*']);
    load_this = find_it(1).name
     data = load(load_this)
     % data  = load_multi_sim(list.name)
    %  data = post_dpf_corrections_frommat(data, load_this);

if ~isfield(rq , 'moyalTest_xyzCorrPA')
    rq.moyalTest_xyzCorrPA = zeros(1, length(rq.CorrS2Chi2_dof));
end


for i = 1:length(rq.CorrS2Chi2_dof)
    if rq.CorrS2Chi2_dof(i)> 0 % if there was a track computed
        
        pulsetypes = data.pulse_classification(: , i);
        x_corr_cm = data.x_corrected(: , i);
        xyz_corr_pulse_area_phe = data.xyz_corrected_pulse_area_all_phe(: , i);
        pulse_filter = (pulsetypes == 2 | pulsetypes == 4  | data.s2_like_class5(: ,i) ==1)  & (x_corr_cm > -99.);
    
        pulse_length = data.pulse_length_samples(:, i );
        pulse_len_per_samp = double(xyz_corr_pulse_area_phe) ./ double(pulse_length);
        
        plps_cut = pulse_len_per_samp > 1;
        pulse_filter = pulse_filter .* plps_cut; %PAT added 0818
        pulse_filter = pulse_filter > 0; % reconvert to logical
        mypulseareapersample = double(xyz_corr_pulse_area_phe(pulse_filter))./ double(pulse_length(pulse_filter));
        
        rq.moyalTest_xyzCorrPA(i) = moyal_test(double(mypulseareapersample));
        
    end
end

    C= strsplit(current_file, '/');
    endString = C{end};
    endString=endString(1:19);
    outName = [ outputpath endString];
    save(outName , 'rq');

end





