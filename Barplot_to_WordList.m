%We want to format our data for the topic model
%need each row to be an observation in the format 
%Timestamp, 0, 0, 0, 1, 2, 2, 4 etc. 
%where the integers are the indeces of the hit repeated as many times as
%there were hits at that observation

%fist input is barplot results, each row corresponding to a taxon, each column
%an observation, and entries are frequency of hits
%second input is date in datetime format
%output Wordlist first column is timestamp in days since first observation
%word_dict relates the word tag to the taxonomic assignment

function [Wordlist, word_dict] = Barplot_to_WordList(Barplot_matrix, dates, index)

bm = Barplot_matrix; %for simplicity
dates = datenum(dates); 

%remove taxonomic assigments that don't have any hits at any observations
rmind = find(sum(bm')==0);
bm(rmind, :) = []; 
word_dict = index; 
word_dict(rmind) = []; 

B = size(bm,1); %-2 %decide if including Eukaryotes and Unassigned, then B = size(bm,1). 

%get dimensions for output
samplesize = sum(bm); 
samplesize = max(samplesize); 
%If bm is already resampled this should be one number. If not, some rows of
%wordlist will have nans at the end. 
num_obs = size(bm, 2); 

%initialize output
Wordlist = nan(num_obs, samplesize+1); %one column for timestamp

T0 = min(dates); 

for obs = 1:num_obs %go through each sample
    %first put in timestamps
    T = dates(obs) - T0; 
    Wordlist(obs, 1) = T;
    
    lastbin = 1; %this is to keep track of how much of row is filled
    for i = 1:B %go through each taxon assignment
        x = bm(i, obs); %x is number of hits
       
        if x > 0
        %insert a vector of index number repeated x times 
        Wordlist(obs, lastbin+1:lastbin+x) = ones(1, x).*(i-1); 
        end
        
        lastbin = lastbin+x; 
        
    end
    
end



end
