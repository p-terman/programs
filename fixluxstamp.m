
inputpath_dp = '/media/paul/TOSHIBA/100_pulse/';
inputpath_rq = '/home/paul/LUXdata/lip_filter_new/live/';
for i = 1 : length(files)   
    current_file = files(i).name
    [dp, lt] = load_multi(current_file , inputpath_dp);
    C = current_file(1:19);
    current_file = [inputpath_rq C '.mat'];
    load(current_file);
    rq.
    for j = 1: length(rq.luxstamp)
        if rq.CorrS2Chi2_dof(j) > 0 
            rq.luxstamp(j) = dp.luxstamp_samples(j);
            
        end
    end
    save_name= ['/home/paul/LUXdata/lip_filter_new/live/with_luxstamp/' C '.mat'];
    save( save_name , 'rq');
end
    