     
for i = 1:length(files)
    current_file = files(i).name
    load(files(i).name);
    current_file = current_file(1:end-4);
    a = strcat(current_file, '_f000000001.rq.mat');
    data = load(['/home/paul/LUXdata/' a]);
    
    
    for j=1:length(rq.total_Corr_keV)
        pulsetypes = data.pulse_classification(: , j);
        pulse_area_phe = data.pulse_area_phe(:,j);
        pulse_length = data.pulse_length_samples(:, j);
        x_corr_cm = data.x_corrected(: , j);
        y_corr_cm = data.y_corrected(: , j);

        drift_us = 10. * 0.001 * double(data.z_drift_samples(: , j));
        drift_to_cm = 1.51/10.;% #cm/us
        z_cm = 50. - drift_us * drift_to_cm  ;
        pulse_filter = (pulsetypes == 2 | pulsetypes == 4 | data.s2_like_class5(: , j) == 1 )  & (x_corr_cm > -99.);

        pulse_len_per_samp = double(pulse_area_phe) ./ double(pulse_length);
        plps_cut = pulse_len_per_samp > 1;
        pulse_filter = pulse_filter .* plps_cut; 
        pulse_filter = pulse_filter > 0; % reconvert to logical
        
        if sum(pulse_filter) > 4 % was >2

            xarr = x_corr_cm(pulse_filter);
            yarr = y_corr_cm(pulse_filter);
            zarr = z_cm(pulse_filter);
            
            data_vec = [xarr(end)-xarr(1), yarr(end)-yarr(1), zarr(end)-zarr(1)];
            data_dir = data_vec ./ norm(data_vec);
            params_opt = rq.CorrS2TrackFit_params (:,j); 
            track_dir = params_opt(4:6); %want the slopes of the parameterized lines
            track_dir = track_dir./norm(track_dir); 
            rq.CorrS2TrackTheta(j) = acos(abs(track_dir(3)))*180./pi;
            rq.CorrS2TrackPhi(j) = atan2(track_dir(2), track_dir(1))*180./pi;
            rq.CorrS2DataToFitAng(j) = acos(dot(data_dir, track_dir)) * 180. / pi;
        end
    end
    save(files(i).name , 'rq');
    
end
