function Run_LIP_filter_on_data(list)
for i = 1:length(list)
    current_file =   string(list(i).name);
    inputpath = '/media/paul/TOSHIBA/matfiles/Original_dpf/';
    outputpath = '/media/paul/TOSHIBA/matfiles/Original_dpf/filtered' ;%lip_filter_new/live/';
    [dp, lt_hour] = load_multi (current_file, inputpath);
    dp = post_dpf_corrections_frommat(dp, current_file);
    rq = LIP_filter_new(dp);
    rq.lt_hour = lt_hour;
    C= strsplit(current_file, '/');
    endString = C{end};
    endString=endString(1:19);
    outName = [ outputpath endString];
    save(outName , 'rq');
    rejoice = ['Done with ' endString '!']
end
    party = ['Done with all live data!']
end
