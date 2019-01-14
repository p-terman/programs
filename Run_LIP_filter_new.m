
inputpath = '/home/paul/matfiles/dp_plus_chop/';
outputpath = '/home/paul/matfiles/dp_plus_chop/filtered/after_moyal/';
%files = dir([inputpath 'luxsm*merged.rq.mat']);
%this is for sim - don't use on data
for k = 1:length(files)
    current_file = [inputpath files(k).name];
    dp = load(current_file);
    %dp = load_multi_sim(current_file);
    rq = LIP_filter_new(dp);
    C= strsplit(current_file, '/');
    endString = C{end};
    endString=endString(1:19);
    outName = [ outputpath endString];
    save(outName , 'rq');
    rejoice = ['Done with ' endString '!']

end

    