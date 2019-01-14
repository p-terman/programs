addpath(genpath(['/media/paul/TOSHIBA/matfiles/live_filtered/after_moyal/']));
      inputpath = '/media/paul/TOSHIBA/100_pulse';
    outputpath = '/media/paul/TOSHIBA/matfiles/live_filtered/refilter/' 
list = dir(['/media/paul/TOSHIBA/matfiles/live_filtered/after_moyal/' 'lux*.mat']);
for pos=1: length(list)
    current_file=list(pos).name;
    load(current_file);
    rq.PAbeforeS1 = zeros(10, length(rq.nS1));
    current_file = current_file(1:19);
    [data, ~] = load_multi_scratch (current_file, inputpath);
  %  data = post_dpf_corrections_frommat(data, current_file);
%no corrections for s1 pulses. 
    data = which_class5_type(data);

     for i = 1: length(rq.nS1)
     n = rq.nS1(i);
     pulsetypes = data.pulse_classification(: , i);
     pulse_filter = pulsetypes == 1 | data.s1_like_class5(:,i) == 1;
     pulse_area_phe = data.pulse_area_phe(:,i);
    if n > 0
        
        %#get the cumulative area before the first 10 S1s
        %#(maybe useful for cutting e-burps)
        pulse_indices_S1s = find(pulse_filter ==1 | data.s1_like_class5(:,i) ==1);
        for ii = 1: length(pulse_indices_S1s)
           ind = pulse_indices_S1s(ii);
            if ii > 10
                break;
            end
            if ind == 1
                rq.PAbeforeS1(ii, i) = 0; %this is for the event - theoretically already zero...
            else
                if ii==1
                    rq.PAbeforeS1(ii , i) = sum(pulse_area_phe(1:(ind-1)));
                else 
                    rq.PAbeforeS1(ii,i) = sum(pulse_area_phe(pulse_indices_S1s(ii - 1) + 1:(ind-1))); % just the difference between this and the last
            
                end                
            end
        end
    end
     end
     save([outputpath current_file '.mat'] , 'rq');
end