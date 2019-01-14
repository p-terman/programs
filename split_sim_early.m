function n_divisions = split_sim_early(filename_evt,data_path_evt,filename_rq,data_path_rq, num_in_each) %,data_processing_settings_path,lug_iqs_xml_file)
% this function replaces 'make_smaller_sim
%this will break up the events of a large sim file into smaller ones. 
%for use with sims NOT LIVE - will end up duplicating livetime \
%- future mods could fix
dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
rq_list = fieldnames(dp);

%num_in_each = 50; % number of events in each new file
n_divisions = ceil(length(dp.event_number)/ num_in_each);

livetime = dp.admin.livetime; %NOTE: this will end up duplicating livetime - doesn't matter for sim
                                %Don't use for live data

 %{
 no cvt yet.
%first do the cvt files. We must do cvt first, because the dp file being
%saved is dependent upon 'settings'. 
[cvt_struct, settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));                                        
cvt_list = fieldnames(cvt_struct);
for k = 1:n_divisions %go through each group for new files
    clear new_cvt        
    
    if length(cvt_struct) -  (k-1)*num_in_each == 0
        break % this should end the loop
    end
    
    if   length(cvt_struct) -  k*num_in_each < 0
        ending = length(cvt_struct) - (k-1) * num_in_each;
    else
        ending = num_in_each;
    end
    
    for i = 1:length(cvt_list)
        current_cvt = cvt_list{i};
        
        for j = 1:ending % go through each event
            
            if j==1
                if strcmp(current_cvt, 'ch_map') || strcmp(current_cvt, 'filename_prefix') ||...
                    strcmp(current_cvt, 'posttrigger') || strcmp(current_cvt, 'pretrigger')
            
                    new_cvt(j).(current_cvt) = cvt_struct(:,1).(current_cvt);
                    
                    
                else
                    new_cvt(j).(current_cvt) = cvt_struct(:,(k-1)*(num_in_each)+j).(current_cvt); 
                end
            elseif j == num_in_each
                if  strcmp(current_cvt, 'filename_prefix') ||...
                    strcmp(current_cvt, 'posttrigger') || strcmp(current_cvt, 'pretrigger')
                    
                    new_cvt(j).(current_cvt) = cvt_struct(:,end).(current_cvt);
                else
                    new_cvt(j).(current_cvt) = cvt_struct(:,(k-1)*(num_in_each)+j).(current_cvt); 
                end    
            else % j not 1 or end
            new_cvt(j).(current_cvt) = cvt_struct(:,(k-1)*(num_in_each)+j).(current_cvt); 
            end
        end
    end % cvt_list

%save the new file
fprintf('Saving CVT file... ');
    filename_cvt_fullpath_new = [data_path_evt filesep filename_evt];
    filename_cvt_fullpath_new = [filename_cvt_fullpath_new(1:end-6) num2str(k) '0.cvt'];
    status = LUXCVTWriter_framework( new_cvt, settings, livetime, filename_cvt_fullpath_new );
end

          
clear cvt_struct
clear new_cvt

%}
                                                           
%now do the dp files                                
for k = 1:n_divisions %go through each group for new files
    clear dp_new
     if length(dp.event_number) -  (k-1)*num_in_each == 0
         break; % this should end the loop
     end    
     
     if   length(dp.event_number) -  k*num_in_each < 0
        ending = length(dp.event_number) - (k-1)*num_in_each;
     else
        ending = num_in_each;
     end
        
    for ii=1:length(rq_list)
        current_rq = rq_list{ii};
        if strcmp(current_rq, 'admin')
            dp_new.(current_rq) = dp.(current_rq);
        elseif strcmp(current_rq, 'source_filename')
            dp_new.(current_rq)( 1:ending , 1:39) ...
                = dp.(current_rq)((k-1)*num_in_each + 1 : (k-1)*num_in_each + ending , 1:39)
         %  = dp.(current_rq)((k-1)*num_in_each + 1 : (k-1)*num_in_each + ending, (k-1)*num_in_each + 1 : (k-1)*num_in_each + ending);
            
        else
            dp_new.(current_rq)(: , 1:ending) ...
                = dp.(current_rq)(: , (k-1)*num_in_each + 1 : (k-1)*num_in_each + ending);
        end
    
    end
    % save the new file
    
    settings.evt_settings = dp.admin.evt_settings;
    settings.daq_settings = dp.admin.daq_settings;
    settings.filename_prefix = dp.admin.filename_prefix;


    fprintf('Saving RQ file... ');
    filename_rq_new = [filename_rq(1:end-5) num2str(k) '0.rq'];
    status = LUXBinaryRQWriter_framework(settings, dp_new, filename_rq_new, data_path_rq, dp_new.event_number, livetime); 
    
end %making the new files

clear dp
clear dp_new

event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);
evt_list = fieldnames(event_struct);

for k = 1:n_divisions %go through each group for new files
    clear new_evt        
    if length(event_struct) - (k-1)*num_in_each == 0
        break % this should end the loop
    end
    
     if   length(event_struct) -  k*num_in_each < 0
            ending = length(event_struct)  - (k-1) * num_in_each;
     else
        ending = num_in_each;
     end
    
    for i = 1:length(evt_list)
        current_evt = evt_list{i};
      
        for j = 1:ending % go through each event
            if j==1
                if strcmp(current_evt, 'ch_map') || strcmp(current_evt, 'filename_prefix') ||...
                    strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')
                    
                    new_evt(j).(current_evt) = event_struct(:,1).(current_evt);
                    
                elseif strcmp(current_evt, 'filename')
                    filename_evt_new = [filename_evt(1:end-6) num2str(k) '0.evt'];
                    new_evt(j).(current_evt) = filename_evt_new;
                else
                    new_evt(j).(current_evt) = event_struct(:,(k-1)*(num_in_each)+j).(current_evt); 
                end
            elseif j == num_in_each
                if  strcmp(current_evt, 'filename_prefix') ||...
                    strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')
                    
                    new_evt(j).(current_evt) = event_struct(:,end).(current_evt);
                else
                    new_evt(j).(current_evt) = event_struct(:,(k-1)*(num_in_each)+j).(current_evt); 
                end
                
            else
                % j not 1 or last 
            new_evt(j).(current_evt) = event_struct(:,(k-1)*(num_in_each)+j).(current_evt); 
            end
        end
    end % evt_list

%save the new file
fprintf('Saving EVT file... ');
    filename_evt_fullpath_new = [data_path_evt filesep filename_evt_new '.mat'];
    save(filename_evt_fullpath_new, 'new_evt');
end


fprintf('Done!');

end

