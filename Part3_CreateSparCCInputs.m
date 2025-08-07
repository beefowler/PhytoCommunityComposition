% This script is for producing the input csvs needed for SparCC analysis

% Primarily, cutting taxonomic assingemnts at Species level with 80%
% confidence, sorting chronologically, and adding a column for whether we%
% consider each taxa a small phytoplankton or a small choroophyte 


%% Load and organize initial data 

%sequence counts 
   load('speciestable.mat') 
%metadata in same order as sample columns 
   sorted_meta = readtable("MVCO_metadata_AllSequenceSamples.csv");


%%  cut down based on confidence at species level 
temp = seqtab_noreps(seqtab_noreps.Species_boot>80, :);
speciestable = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
speciestable.total = sum(speciestable{:,10:134},2);
speciestable = sortrows(speciestable, 'total', 'descend');

% sort chronologically 
[a, b]= sort(meta_noreps.sampledates); 
temp = 10:134;
chrono_speciestable = speciestable(:, [1:8 temp(b)]); %now chronological 



%%  Now add column for taxa we consider picos and chlorophytes

chrono_speciestable.Picos = ones(height(chrono_speciestable), 1); 

% and a column for pico chlorophytes 

chrono_speciestable.Picos = ones(height(chrono_speciestable), 1); 

    %(This is copied and modified from Part2 cutting function, so make sure
    %any changes there are implemented here) 
    
chrono_speciestable.Picos(~strcmp(speciestable.Kingdom, 'Eukaryota'), :) = 0;
chrono_speciestable.Picos(strcmp(speciestable.Supergroup, 'Opisthokonta'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Supergroup, 'Amoebozoa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Supergroup, 'Apusozoa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Apicomplexa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Cercozoa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Pseudofungi'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Foraminifera'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Ciliophora'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Discoba'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Katablepharidophyta'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Metamonada'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Opalozoa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Perkinsea'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Picozoa'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Rhodophyta'), :) = 0;
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Sagenista'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Streptophyta'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Telonemia'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Colpodellidea'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Diplonemea'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Embryophyceae'), :) = 0;
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Heterolobosea'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Kinetoplastea'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Labyrinthulomycetes'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Syndiniales'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Ellobiopsidae'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Family, 'Spironemidae'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Family, 'Rhodelphidae'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Filipodium'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Rapaza'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Phyllosiphon'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Ulvophyceae'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Pellita'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Centroheliozoa'), :) = 0;
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Radiolaria'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Euglenida'), :) = 0;
chrono_speciestable.Picos(strcmp(speciestable.Family, 'Distigmidae'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Division, 'Dinoflagellata'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Ansanella_granifera'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Amphidinium_carterae'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Biecheleria_cincta'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Biecheleriopsis_adriatica'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Gymnodinium_smaydae'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Heterocapsa_rotunda'), :) = 1;
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Heterocapsa_circularisquama'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Karlodinium_veneficum'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Lepidodinium_viride'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Paragymnodinium'), :) = 1;
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Prorocentrum_cordatum'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Genus, 'Yihiella'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Class, 'Bacillariophyta'), :) = 0; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Arcocellulus_cornucervis'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Chaetoceros_throndsenii'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Minidiscus_comicus'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Minidiscus_trioculatus'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Minutocellus_polymorphus'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Skeletonema_grethae'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Skeletonema_japonicum'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Skeletonema_marinoi'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Skeletonema_menzellii'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Skeletonema_pseudocostatum'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Thalassiosira_mala'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Thalassiosira_oceanica'), :) = 1; 
chrono_speciestable.Picos(strcmp(speciestable.Species, 'Thalassiosira_pseudonana'), :) = 1; 



chrono_speciestable.PicosChlorophytes = chrono_speciestable.Picos;
chrono_speciestable.PicosChlorophytes(~strcmp(speciestable.Division, 'Chlorophyta'), :) = 0;

taxastring = chrono_speciestable.Kingdom; 
for i = 1:height(chrono_speciestable)
    string = taxastring{i};
    string = [string ',' chrono_speciestable.Supergroup{i} ',' chrono_speciestable.Division{i} ',' chrono_speciestable.Class{i} ',' chrono_speciestable.Order{i} ',' chrono_speciestable.Family{i}  ',' chrono_speciestable.Genus{i}  ',' chrono_speciestable.Species{i}];
    taxastring{i} = string; 
end
clear string
chrono_speciestable.taxastring = taxastring; % add this string to table 


% save 
writetable(chrono_speciestable, 'chrono_species_gt80.csv')