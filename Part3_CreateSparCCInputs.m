% This script is for producing the input csvs needed for SparCC analysis

% Primarily, cutting taxonomic assingemnts at Species level with 80%
% confidence, sorting chronologically, and adding a column for whether we%
% consider each taxa a small phytoplankton or a small choroophyte 


%% Load data 

%processed sequencing table 
load("seqtab_noreps.mat")


%%  cut down based on confidence at species level, keeping track of which are small phyto and chlorophytes

is_chlorophyta = strcmp(small_phyto_noreps.Division, 'Chlorophyta'); 
small_phyto_noreps.is_chlorophyta = is_chlorophyta; 

temp = small_phyto_noreps(small_phyto_noreps.Species_boot>80, :);
all_Species = grpstats(small_phyto_noreps, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', small_phyto_noreps.Properties.VariableNames([19:143 145 147]));

%Convert boolean vectors back to 1 or 0 
all_Species.sum_small_phytoplankton = all_Species.sum_small_phytoplankton>0;
all_Species.sum_is_chlorophyta = all_Species.sum_is_chlorophyta>0;

all_Species.total = sum(all_Species{:,10:134},2);
all_Species = sortrows(all_Species, 'total', 'descend');

% sort chronologically 
[a, b]= sort(meta_noreps.sampledates); 
temp = 10:134;
chrono_speciestable = all_Species(:, [1:8 temp(b) 135 136]); %now chronological

% save 
writetable(chrono_speciestable, 'chrono_species_gt80.csv')