%% We want to make two Wordlists for input into topic modeling software for small phytoplankton and chlorophyte species respectively

% Start as before by loading data


%% Load data 

%processed sequencing table 
load("seqtab_noreps.mat")


%%  cut down based on confidence at species level 
temp = small_phyto_noreps(small_phyto_noreps.Species_boot>80, :);
small_phyto_Species = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
small_phyto_Species.total = sum(small_phyto_Species{:,10:134},2);
small_phyto_Species = sortrows(small_phyto_Species, 'total', 'descend');

temp = cut_to_chlorophytes(temp); %further cut to only cholorophyta
chlorophyta_Species = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
chlorophyta_Species.total = sum(chlorophyta_Species{:,10:134},2);
chlorophyta_Species = sortrows(chlorophyta_Species, 'total', 'descend');

%%  sort chronologically 
[a, b]= sort(meta_noreps.sampledates); 
temp = 10:134;
chrono_chloro_species = chlorophyta_Species(:, [1:9 temp(b)]); %now chronological 
chrono_pico_species = small_phyto_Species(:, [1:9 temp(b)]); %now chronological 
chrono_meta_noreps = meta_noreps(b, :); 
dates = unique(meta_noreps.sampledates);

%% Convert to WordLists & Save as CSV 

[Wordlist, word_dict] = Barplot_to_WordList(table2array(chrono_chloro_species(:,10:end)), dates, chrono_chloro_species.Properties.RowNames);

wordlist_cell = num2cell(Wordlist); 
mask = cellfun(@(C) all(isnan(C)), wordlist_cell); 
wordlist_cell(mask) = {''};

save('Chlorophyta_wordlist.mat', 'Wordlist', 'word_dict')

writecell(wordlist_cell, 'new_4Batch_chlorophyta.csv', 'FileType', 'text');



[Wordlist, word_dict] = Barplot_to_WordList(table2array(chrono_pico_species(:,10:end)), dates, chrono_pico_species.Properties.RowNames);

wordlist_cell = num2cell(Wordlist); 
mask = cellfun(@(C) all(isnan(C)), wordlist_cell); 
wordlist_cell(mask) = {''};

save('SmallPhyto_wordlist.mat', 'Wordlist', 'word_dict')
writecell(wordlist_cell, 'new_4Batch_SmallPhytos.csv', 'FileType', 'text');



%% Here is the function where we can change which sepcies we want included in small phytoplankton  




function cutable = cut_to_chlorophytes(speciestable)

cutable = speciestable; 
cutable(~strcmp(cutable.Division, 'Chlorophyta'), :) = [];

end