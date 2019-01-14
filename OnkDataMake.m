db=readtable('/home/paul/Downloads/OnkelosGenesis-he-OnkelosGenesis.xlsx')
%this will make db with first line wrong
% convert to cell using db = table2cell(db);
%then copy paste the first line in 
%then transfer to dababase by hand because we don't have the special
%package

new_count=0;
for i= 1:size(db,1)
    verse= db(i,2);
    verse=verse{1}; %convert from cell to string
    
    words=strsplit(verse);
    len=length(words);
    
    info=db(i,1);
    info=info{1};
    
    parse_info = strsplit(info);
    book=parse_info(2);
    citation=parse_info(3);
    citation=citation{1};
    citation=strsplit(citation, ':');
    chapter=citation(1);
    verse_num=citation(2);
    wordnr=0;
    for j=1:len
        wordnr=wordnr+1;
        new_count=new_count+1;  
     
        new_Array(new_count, 1) = {new_count};
        new_Array(new_count, 2) = book;
        new_Array(new_count, 3) = chapter;
        new_Array(new_count, 5) = {wordnr};
        new_Array(new_count, 4) = verse_num;   
        new_Array(new_count, 6) = words(j);
    end
end

%{ 
Might be useful for removing/combining words
for i=2:17
value = new_Array{i,4};
value=value-1;
new_Array(i,4)={value};
end
will need similar program for (i,1) as well but going throughout the
db
%}