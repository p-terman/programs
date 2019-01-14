function  [total_struct , n_total] = Load_filtered_LIPs_improved

clear total_struct; 
files = dir("lux*.mat") % or we could force a speific list as input instead. 

%clear data_struct;

n_total.total = 0;
n_total.hasfit = 0;
n_total.has5S2 = 0;

load(files(1).name);
rq_list = fieldnames(rq); 
    
for j = 1:length(rq_list) 
    
    current_rq = rq_list{j};
    if strcmp(current_rq, 'lt_hour')
        total_struct.lt_hour = 0;
    else
        total_struct.(current_rq) = [];
    end
end

for i = 1:length(files)
 %   struct_name = files(i).name(11: end-4);

    load(files(i).name) % this should load rq 

    for j = 1:length(rq_list)
        current_rq = rq_list{j};
        j
        if strcmp(current_rq, 'lt_hour')
            total_struct.lt_hour = total_struct.lt_hour + rq.lt_hour;
        else
            total_struct.(current_rq) = [total_struct.(current_rq) , rq.(current_rq)(rq.CorrS2Chi2_dof > 0) ];
        end    
    end
  
    n_total.total = n_total.total + length(rq.CorrS2Chi2_dof);
    n_total.hasfit = n_total.hasfit + length(rq.CorrS2Chi2_dof(rq.CorrS2Chi2_dof>0) == 1); 
    n_total.has5S2 = n_total.has5S2 + length(rq.nCorrS2PulsesBeforeCut(rq.nCorrS2PulsesBeforeCut > 4) == 1);
end

end
    
    
        