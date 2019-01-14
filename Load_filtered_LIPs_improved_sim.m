function  [total_struct , n_total] = Load_filtered_LIPs_improved_sim

clear total_struct n_total;
%files = dir("lux*.mat") % or we could force a speific list as input instead. 
%charge_list = [ 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.15 0.20 0.25 0.30];
charge_list = [0.05 0.10]
%clear data_struct;

n_total.total.total = 0;
n_total.hasfit.total = 0;
n_total.has5S2.total = 0;
n_total.hasEdep.total = 0;
n_total.hasfit.all = zeros(1,length(charge_list));
n_total.has5S2.all = zeros(1,length(charge_list));
n_total.hasEdep.all = zeros(1,length(charge_list));
n_total.hasEvt.all = zeros(1,length(charge_list));

for k = 1:length(charge_list)

    current_file = string(charge_list(k));
    C= strsplit(current_file, '.');
    if length(char(current_file)) == 3
          get_these = strcat("*" , C(2) , "00.mat");
    else
          get_these = strcat("*" , C(2) , "0.mat");
    end
    
    files = dir(get_these);
    load(files(1).name);
    rq_list = fieldnames(rq); 
    
    name_for_charge = strcat('s', C(2));
    n_total.total.(name_for_charge) = 0;    
    n_total.hasfit.(name_for_charge) = 0;
    n_total.has5S2.(name_for_charge) = 0;
    n_total.hasEdep.(name_for_charge) = 0;
    n_total.hasEvt.(name_for_charge) = 0;
     if k == 1
        for j = 1:length(rq_list) 
            current_rq = rq_list{j};
            if strcmp(current_rq, 'lt_hour')
                total_struct.lt_hour = 0;
            else
                total_struct.(current_rq) = [];
            end
        end
        total_struct.charge = [];
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
    
    total_struct.charge(end+1: end+length(rq.CorrS2Chi2_dof(rq.CorrS2Chi2_dof>0))) =  charge_list(k);
    n_total.total.total = n_total.total.total + length(rq.CorrS2Chi2_dof);
    n_total.hasfit.total = n_total.hasfit.total + length(rq.CorrS2Chi2_dof(rq.CorrS2Chi2_dof>0) == 1); 
    n_total.has5S2.total = n_total.has5S2.total + length(rq.nCorrS2PulsesBeforeCut(rq.nCorrS2PulsesBeforeCut > 4) == 1);
    n_total.hasEdep.total = n_total.hasEdep.total + length(rq.raw_total_keV(rq.raw_total_keV>0));
    
    n_total.hasEvt.(name_for_charge) = n_total.hasEvt.(name_for_charge) +  length(rq.CorrS2Chi2_dof);
    n_total.hasfit.(name_for_charge) = n_total.hasfit.(name_for_charge) +  length(rq.CorrS2Chi2_dof(rq.CorrS2Chi2_dof>0) == 1);
    n_total.has5S2.(name_for_charge) = n_total.has5S2.(name_for_charge) + length(rq.nCorrS2PulsesBeforeCut(rq.nCorrS2PulsesBeforeCut > 4) == 1);
    n_total.hasEdep.(name_for_charge) = n_total.hasEdep.(name_for_charge) + length(rq.raw_total_keV(rq.raw_total_keV>0));

end
    n_total.total.(name_for_charge) = n_total.hasfit.(name_for_charge);
    n_total.hasfit.all(k) = n_total.hasfit.(name_for_charge);
    n_total.has5S2.all(k) = n_total.has5S2.(name_for_charge);
    n_total.hasEdep.all(k) = n_total.hasEdep.(name_for_charge);
    n_total.hasEvt.all(k) = n_total.hasEvt.(name_for_charge);
end
    
end
        