% This script picks up with analysis from the outputs of the R processesing

%% Load and organize initial data 

%sequence counts 
   seqtab = readtable("mergtab4_trans.txt"); 
%metadata in same order as sample columns 
   sorted_meta = readtable("MVCO_metadata_AllSequenceSamples.csv");

%some names got messed up. Fix them here 
flist = find(isnan(seqtab.Kingdom_boot));
for i = flist'
    seqtab(i, 9:15) = seqtab(i, 10:16); 
    temp = seqtab{i,17}; 
    seqtab{i, 16} = str2num(temp{:});
    seqtab(i,17) = seqtab(i,18); 
    seqtab{i,18} = {''};
    seqtab(i,19:end-1) = seqtab(i,20:end); 
end
seqtab(:,end) = [];


%want to make a single string that has complete taxa names
taxastring = seqtab.Kingdom; 
for i = 1:height(seqtab)
    string = taxastring{i};
    string = [string ',' seqtab.Supergroup{i} ',' seqtab.Division{i} ',' seqtab.Class{i} ',' seqtab.Order{i} ',' seqtab.Family{i}  ',' seqtab.Genus{i}  ',' seqtab.Species{i}];
    taxastring{i} = string; 
end
clear string
seqtab.taxastring = taxastring; % add this string to seqtab table 


% Remove replicates from sequences and metadata table 
indtoinclude = zeros(height(sorted_meta), 1); 
indtoinclude(~strcmp(sorted_meta.size_fraction, '8_um')) = 1;
dates = unique(sorted_meta.sampledates);
for i = 1:length(dates)
    a = find(indtoinclude & sorted_meta.sampledates == dates(i));
    if length(a) > 1 
        disp(a)
        indtoinclude(a(1)) = 0; 
    end
end

meta_noreps = sorted_meta(logical(indtoinclude), :); 
seqtab_noreps = seqtab(:, [1:18 18+find(indtoinclude') 153]);

clear seqtab taxastring
clear sorted_meta 
clear indtoinclude i a ans flist temp

%% Cut ASV list to small eukaryote phytoplankton only 

small_phyto_vec = ones(height(seqtab_noreps),1); 
small_phyto_reason = strings(height(seqtab_noreps),1); 

% First cut prokaryotes  
small_phyto_vec(~strcmp(seqtab_noreps.Kingdom, 'Eukaryota'), :) = 0; 
small_phyto_reason(~strcmp(seqtab_noreps.Kingdom, 'Eukaryota'), :) = "Not Eukaryote"; 


% Then cut things that aren't phytoplankton (heterotrophs, land plants, parasites etc.)
small_phyto_vec(strcmp(seqtab_noreps.Supergroup, 'Opisthokonta'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Supergroup, 'Opisthokonta'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Supergroup, 'Amoebozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Supergroup, 'Amoebozoa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Supergroup, 'Apusozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Supergroup, 'Apusozoa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Apicomplexa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Apicomplexa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Cercozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Cercozoa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Pseudofungi'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Pseudofungi'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Foraminifera'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Foraminifera'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Ciliophora'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Ciliophora'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Discoba'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Discoba'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Katablepharidophyta'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Katablepharidophyta'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Metamonada'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Metamonada'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Opalozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Opalozoa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Perkinsea'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Perkinsea'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Picozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Picozoa'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Rhodophyta'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Rhodophyta'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Sagenista'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Sagenista'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Streptophyta'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Streptophyta'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Telonemia'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Telonemia'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Colpodellidea'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Colpodellidea'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Diplonemea'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Diplonemea'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Embryophyceae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Embryophyceae'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Heterolobosea'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Heterolobosea'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Kinetoplastea'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Kinetoplastea'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Labyrinthulomycetes'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Labyrinthulomycetes'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Syndiniales'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Syndiniales'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Ellobiopsidae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Ellobiopsidae'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Family, 'Spironemidae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Family, 'Spironemidae'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Family, 'Rhodelphidae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Family, 'Rhodelphidae'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Filipodium'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Filipodium'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Rapaza'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Rapaza'), :) = 'Not Phytoplankton';

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Phyllosiphon'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Phyllosiphon'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Ulvophyceae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Ulvophyceae'), :) = 'Not Phytoplankton'; 

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Pellita'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Pellita'), :) = 'Not Phytoplankton'; 


% Now cut groups that are too large to be measured by flow cytometery
small_phyto_vec(strcmp(seqtab_noreps.Division, 'Centroheliozoa'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Centroheliozoa'), :) = 'Not Phytoplankton; Not Target Cell Size'; 

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Radiolaria'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Radiolaria'), :) = 'Not Phytoplankton; Not Target Cell Size'; 

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Euglenida'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Euglenida'), :) = 'Not Phytoplankton; Not Target Cell Size'; 

small_phyto_vec(strcmp(seqtab_noreps.Family, 'Distigmidae'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Family, 'Distigmidae'), :) = 'Not Target Cell Size'; 



% Remove dinoflagellates unless in the following genera with known small cell size

small_phyto_vec(strcmp(seqtab_noreps.Division, 'Dinoflagellata'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Division, 'Dinoflagellata'), :) = 'Not Target Cell Size'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Ansanella_granifera'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Ansanella_granifera'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Amphidinium_carterae'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Amphidinium_carterae'), :) = 'Confirmed Nano-Dinoflagellates';  

%small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Azadinium'), :) = 1; 
%small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Azadinium'), :) = 'Possible Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Biecheleria_cincta'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Biecheleria_cincta'), :) = 'Confirmed Nano-Dinoflagellates';  

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Biecheleriopsis_adriatica'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Biecheleriopsis_adriatica'), :) = 'Confirmed Nano-Dinoflagellates'; 

%small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Effrenium'), :) = 1; 
%small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Effrenium'), :) = 'Possible Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Gymnodinium_smaydae'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Gymnodinium_smaydae'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Heterocapsa_rotunda'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Heterocapsa_rotunda'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Heterocapsa_circularisquama'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Heterocapsa_circularisquama'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Karlodinium_veneficum'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Karlodinium_veneficum'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Lepidodinium_viride'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Lepidodinium_viride'), :) = 'Confirmed Nano-Dinoflagellates'; 

%small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Luciella'), :) = 1; 
%small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Luciella'), :) = 'Possible Nano-Dinoflagellates';  

%small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Nusuttodinium'), :) = 1; 
%small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Nusuttodinium'), :) = 'Possible Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Paragymnodinium'), :) = 1; % all species in our table are confirmed
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Paragymnodinium'), :) = 'Confirmed Nano-Dinoflagellates';  

%small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Pelagodinium'), :) = 1; 
%small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Pelagodinium'), :) = 'Possible Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Prorocentrum_cordatum'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Prorocentrum_cordatum'), :) = 'Confirmed Nano-Dinoflagellates'; 

small_phyto_vec(strcmp(seqtab_noreps.Genus, 'Yihiella'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Genus, 'Yihiella'), :) = 'Confirmed Nano-Dinoflagellates'; 


% Remove diatoms unless in the following genera with known small cell size

small_phyto_vec(strcmp(seqtab_noreps.Class, 'Bacillariophyta'), :) = 0; 
small_phyto_reason(strcmp(seqtab_noreps.Class, 'Bacillariophyta'), :) = 'Not Target Cell Size'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Arcocellulus_cornucervis'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Arcocellulus_cornucervis'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Chaetoceros_throndsenii'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Chaetoceros_throndsenii'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Minidiscus_comicus'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Minidiscus_comicus'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Minidiscus_trioculatus'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Minidiscus_trioculatus'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Minutocellus_polymorphus'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Minutocellus_polymorphus'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Skeletonema_grethae'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Skeletonema_grethae'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Skeletonema_japonicum'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Skeletonema_japonicum'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Skeletonema_marinoi'), :) = 1; 
small_phyto_reason( strcmp(seqtab_noreps.Species, 'Skeletonema_marinoi'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Skeletonema_menzellii'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Skeletonema_menzellii'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Skeletonema_pseudocostatum'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Skeletonema_pseudocostatum'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Thalassiosira_mala'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Thalassiosira_mala'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Thalassiosira_oceanica'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Thalassiosira_oceanica'), :) = 'Confirmed Nano-Diatom'; 

small_phyto_vec(strcmp(seqtab_noreps.Species, 'Thalassiosira_pseudonana'), :) = 1; 
small_phyto_reason(strcmp(seqtab_noreps.Species, 'Thalassiosira_pseudonana'), :) = 'Confirmed Nano-Diatom'; 


seqtab_noreps.small_phytoplankton = small_phyto_vec; 
seqtab_noreps.small_phytoplankton_reason = small_phyto_reason; 

writetable(seqtab_noreps, 'seqtab_noreps.csv')

small_phyto_noreps = seqtab_noreps(logical(seqtab_noreps.small_phytoplankton), :);

save('seqtab_noreps.mat', 'seqtab_noreps', 'meta_noreps', 'small_phyto_noreps')


small_phyto_Species = grpstats(small_phyto_noreps, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', small_phyto_noreps.Properties.VariableNames(19:143));


writetable(small_phyto_Species, 'Taxa_Included_in_SmallPhytoplankton.csv')

clear small_phyto_vec small_phyto_reason 

%% Report numbers of ASVs and look at small phyto grouped at Division 

disp('ASVs in small phyto:')
disp(height(small_phyto_noreps))

temp = cut_to_chlorophytes(small_phyto_noreps);
disp('ASVs in Chlorophyta:')
disp(height(temp)) 

% cut down based on confidence at division level 
temp = small_phyto_noreps(small_phyto_noreps.Division_boot>80, :);
small_phyto_Division = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));

disp("ASV Proportions: ")
disp(small_phyto_Division.Division)
disp(small_phyto_Division.GroupCount./sum(small_phyto_Division.GroupCount))

disp("Read Proportions: ")
small_phyto_Division.total = sum(small_phyto_Division{:,5:129},2);
disp(small_phyto_Division.total./sum(small_phyto_Division.total))

%

figure 
subplot(2,2,1)
labels = small_phyto_Division.Division; 
explode = zeros(length(labels),1);
explode(1) = 1; 
pie(small_phyto_Division.total, explode, labels)
title('Reads')
subplot(2,2,2)
pie(small_phyto_Division.GroupCount, explode, labels)
title("ASV's")

t = [0.5529    0.8275    0.7804;
    1.0000    1.0000    0.7020;
    0.7451    0.7294    0.8549;
    0.9843    0.5020    0.4471;
    0.5020    0.6941    0.8275;
    0.9922    0.7059    0.3843;
    0.7020    0.8706    0.4118;
    0.9882    0.8039    0.8980;
    0.8510    0.8510    0.8510;
    0.7373    0.5020    0.7412;
    0.8000    0.9216    0.7725];
colormap(t)


% Same calculations as above, but look at chlorophyta grouped at Class 

% cut down to based on confidence at class level 
temp = small_phyto_noreps(small_phyto_noreps.Class_boot>80, :);
temp = cut_to_chlorophytes(temp); %further cut to only cholorophyta
chlorophyta_Class = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));

disp("ASV Proportions: ")
disp(chlorophyta_Class.Class)
disp(chlorophyta_Class.GroupCount./sum(chlorophyta_Class.GroupCount))

disp("Read Proportions: ")
chlorophyta_Class.total = sum(chlorophyta_Class{:,6:129},2);
disp(chlorophyta_Class.total./sum(chlorophyta_Class.total))

subplot(2,2,3)
labels = chlorophyta_Class.Class; 
explode = zeros(length(labels),1);
pie(chlorophyta_Class.total, explode, labels)
title('Reads')
subplot(2,2,4)
pie(chlorophyta_Class.GroupCount, explode, labels)
title("ASV's")

    colormap(t)

    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    clear temp explode labels 

%%  Calculate Shannon Diversity Indeces at Genus level 

%start with genus table 
temp = small_phyto_noreps(small_phyto_noreps.Genus_boot>80, :);
small_phyto_Genus = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));


%calculate H for small phytoplankton 
mat = table2array(small_phyto_Genus(:, 9:133)); 
shanon_smallphyto = nan(1,width(mat));
for i = 1:width(mat)
    sample = mat(:,i)+1; 
    N = sum(sample); 
    H = 0; 
    for j = 1:length(sample)
       H = H + sample(j)/N.*log(sample(j)/N);
    end
    shanon_smallphyto(i) = -H; 
end

disp('IQR Small phytoplankton shannon diversity: ')
disp([quantile(shanon_smallphyto, .25) quantile(shanon_smallphyto, .75)])

%

%calculate H for chlorophyta
temp = cut_to_chlorophytes(temp); %cut to small phytoplankton 
chlorophyta_Genus = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));

mat = table2array(chlorophyta_Genus(:, 9:133)); 

shanon_chlorophytes = nan(1,width(mat));
for i = 1:width(mat)
    sample = mat(:,i)+1; 
    N = sum(sample); 
    H = 0; 
    for j = 1:length(sample)
       H = H + sample(j)/N.*log(sample(j)/N);
    end
    shanon_chlorophytes(i) = -H; 
end

disp('IQR chlorophyta shannon diversity: ')
disp([quantile(shanon_chlorophytes, .25) quantile(shanon_chlorophytes, .75)])

clear N H mat temp i j sample t 

%% Supplementary Figure S2

% Barplot prep (get colors ready and sort things in the order we want) 

map = brewermap(12, 'paired');
temp = map(1,:);
map([1 4 5 8 9 12], :) = [];
map1 = brewermap(12, 'Set3');
map(4, :) = map1(end,:);
map(6, :) = temp;


clear temp map1

%
% Get name of dominant genus in every sample 
top_genera = string; 
for i = 9:133
    [sample, orderr] = sort(table2array(small_phyto_Genus(:,i))); 
    top_genera(i) = small_phyto_Genus.Genus(orderr(end)); 
end
top_genera(1:8) = []; %remove empty first entries
[G, l] = findgroups(top_genera); 
[a] = hist(G, 1:max(G));
[a2, ind2] = sort(a, 'descend'); %a2 is now the number of samples dominated by each group in order l(ind2)

boxgroups = [G == ind2(1); G == ind2(2); G == ind2(3); G == ind2(4); G == ind2(5)];
boxgroups = [boxgroups; (sum(boxgroups) == 0)];
%Each row of boxgroups is an idex logical for the 6 boxes we want in the
%supplementary figure 

%get total read numbers
small_phyto_Genus.total = sum(small_phyto_Genus{:,9:133},2);
small_phyto_Genus = sortrows(small_phyto_Genus, 'total', 'descend');

%figure 
% now want a boxplot with diversity for samples with dominant taxa in each of top 6 genera 
subplot(1, 2,1)
cla
order1 = [1 3 4 2 5]; 
for i = 1:5
    hold on 
    g = boxchart(zeros(1, sum(boxgroups(i,:), 2))+i, shanon_smallphyto(boxgroups(i,:)));
    g.JitterOutliers = 'off';
    g.MarkerStyle = '.';
    g.BoxFaceColor = map(order1(i),:);
    g.BoxFaceAlpha = .75;
    g.MarkerColor = 'k';
    g.BoxEdgeColor = 'k';
end
i = 6;
    g = boxchart(zeros(1, sum(boxgroups(i,:), 2))+i, shanon_smallphyto(boxgroups(i,:)));
    g.JitterOutliers = 'off';
    g.MarkerStyle = '.';
    g.BoxFaceColor = 'k';
    g.MarkerColor = 'k';
xticks(1:10)
xticklabels([l(ind2(1:5)) 'Other'])
ylabel('Shannon diversity index')
title('Sample Diversity by Dominant Taxa')
xlim([0 7])

% Check statistical significance 
[h, p] = ttest2(shanon_smallphyto(boxgroups(1,:)), shanon_smallphyto(boxgroups(2,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(1,:)), shanon_smallphyto(boxgroups(3,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(1,:)), shanon_smallphyto(boxgroups(4,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(1,:)), shanon_smallphyto(boxgroups(5,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(2,:)), shanon_smallphyto(boxgroups(3,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(2,:)), shanon_smallphyto(boxgroups(4,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(2,:)), shanon_smallphyto(boxgroups(5,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(3,:)), shanon_smallphyto(boxgroups(4,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(3,:)), shanon_smallphyto(boxgroups(5,:)))
[h, p] = ttest2(shanon_smallphyto(boxgroups(4,:)), shanon_smallphyto(boxgroups(5,:)))


% also want a boxplot with diversity in each month 
subplot(1,2,2)
cla
h = boxchart(month(meta_noreps.sampledates), shanon_smallphyto);
h.BoxFaceColor = 'k';
    h.JitterOutliers = 'off';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';
xlim([0 13])
xticks(1:2:12)
xticklabels({'Jan'; 'Mar'; 'May'; 'Jul'; 'Sep'; 'Nov'})
ylabel('Shannon diversity index')
title('Sample Diversity by Month')

    set(findall(gcf,'-property','FontSize'),'FontSize',16)

    clear g G h i ind2 l p order1 orderr map1 a a2 boxgroups sample top_genera
%% Barplots for Figure 1

% Get name of dominant genus in every sample 
top_genera = string; 
for i = 9:133
    [sample, orderr] = sort(table2array(small_phyto_Genus(:,i))); 
    top_genera(i) = small_phyto_Genus.Genus(orderr(end)); 
end
top_genera(1:8) = []; %remove empty first entries
[G, l] = findgroups(top_genera); 
[a] = hist(G, 1:max(G));
[a2, ind2] = sort(a, 'descend'); %a2 is now the number of samples dominated by each group in order l(ind2)

%get total read numbers
small_phyto_Genus.total = sum(small_phyto_Genus{:,9:133},2);
small_phyto_Genus = sortrows(small_phyto_Genus, 'total', 'descend');

figure 

%First we want a time series of most abundant taxa by read numbers 
newmat = [table2array(small_phyto_Genus(:, 9:133))] ;
newnewmat = [newmat(1:5,:); sum(newmat(7:end,:),1)];

subplot(2,1,1)
alp = bar(meta_noreps.sampledates, ones(1,125), 20, 'stacked');
alp.FaceColor = [.8 .8 .8];
alp.FaceAlpha = 0.6
alp.EdgeColor = 'none'; 

hold on 
a = bar(meta_noreps.sampledates, newnewmat(1:5,:)./sum(newmat), 20, 'stacked');
for i = 1:5
    a(i).FaceColor = map(i,:); 
    a(i).FaceAlpha= (0.85); 
    a(i).EdgeColor = 'none'; 
end
ylabel('Proportion of reads', 'FontSize', 24)
title('Small Phytoplankton')
legend({'Other'; 'Micromonas'; 'Bathycoccus'; 'Picochlorum';'Phaeocystis'; 'Ostreococcus'})
xlim([min(meta_noreps.sampledates) - 70, max(meta_noreps.sampledates) + 70])
small_phyto_Genus = sortrows(small_phyto_Genus, 'total', 'descend');


subplot(2,2,3)
b = bar(6:10, small_phyto_Genus.total(6:10)); 
hold on 
b.FaceColor = [.8 .8 .8]; 
order1 = [1 2 3 4 5 5 5 5 5 5]  %just for coloring 
for i = [1 2 3 4 5] %7 10]
    b = bar(i, small_phyto_Genus.total(i));
    b.FaceColor = map(order1(i),:); 
    b.FaceAlpha= .8;
end

xticks(1:10)
xticklabels(small_phyto_Genus.Genus(1:10))
ylabel('Number of reads', 'FontSize', 24)
title('Read Counts per Genus')

small_phyto_Genus = sortrows(small_phyto_Genus, 'total', 'descend');

subplot(2,2,4)
b = bar(6:10, a2(6:10));
hold on 
b.FaceColor = [.8 .8 .8];  %just for coloring
order1 = [1 3 4 2 5 0 0 0 0 0]; 
for i = [1 2 3 4 5 ]
    b = bar(i, a2(i));
    b.FaceColor = map(order1(i),:); 
    b.FaceAlpha= .8;
end
xticks(1:10)
xticklabels(l(ind2(1:10)))%, 'Interpreter', 'none')
ylim([0 32])
ylabel('Number of samples', 'FontSize', 24)
title('Dominant Genus by Sample')

    set(findall(gcf,'-property','FontSize'),'FontSize',16)

clear a b alp newmat newnewmat l i G a2 order1 orderr sample ind2 top_genera

    %% PCA at Species level assignents - Figure 3

% cut down based on confidence at species level 
temp = small_phyto_noreps(small_phyto_noreps.Species_boot>80, :);

small_phyto_Species = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
small_phyto_Species.total = sum(small_phyto_Species{:,10:134},2);
small_phyto_Species = sortrows(small_phyto_Species, 'total', 'descend');

temp = cut_to_chlorophytes(temp); %further cut to only cholorophyta
chlorophyta_Species = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
chlorophyta_Species.total = sum(chlorophyta_Species{:,10:134},2);
chlorophyta_Species = sortrows(chlorophyta_Species, 'total', 'descend');

%Now use these tables to calculate CLR 
rawvalues = table2array(chlorophyta_Species(:,10:134));
n = height(rawvalues); 
offset = rawvalues+1; 
CLR = log(offset) - (1/n) .* sum(log(offset)); 

[wcoeff_chlorophytes,score,latent,tsquared,explained] = pca(CLR'); %this is the crux of the PCA


figure
subplot(4,2,2)

disp('Variance Explained:')
disp(explained(1:2))

%plot a scatter plot of PCA 
scatter(score(:,1), score(:,2), 35, day(meta_noreps.sampledates, 'dayofyear'), 'filled')
xlabel('Principal component 1 (21.7%)')
ylabel('Principal component 2 (17.0%)')
title('Chlorophyta')
caxis([0 366])
ylim([-13 13])
xlim([-11 12])
box on 

map = brewermap(11,'spectral');
map(3:6,:) = []; 
map = [map(3:end, :); map(1:2, :)];
map = map(end:-1:1, :);
colormap(map)
%to make warm colors in summer cool in winter
map = [map(4:end, :); map(1:3,:)]; 
map(4,:) = [0.94  0.608  0.361]; 
colormap(map)

h = colorbar;
colormap(map);
caxis([0 366])
h.Ticks = [0 91 166 258]; 
h.TickLabels = {'Jan 01'; 'Apr 01'; 'Jun 15'; 'Sep 15'}; 
set( h, 'YDir', 'reverse' );

%Plot a time series 
subplot(4,1,3)
[~, chrono_index] = sort(meta_noreps.sampledates); %get chronological index so we can plot connected dots 
cla
plot(meta_noreps.sampledates(chrono_index), score(chrono_index,2), 'color', [.8 .8 .8], 'linewidth', 2)
hold on 
scatter(meta_noreps.sampledates(chrono_index), score(chrono_index,2), 20, [.8 .8 .8], 'filled')
plot(meta_noreps.sampledates(chrono_index), score(chrono_index,1), 'color', [.2 .2 .2], 'linewidth', 2)
scatter(meta_noreps.sampledates(chrono_index), score(chrono_index,1), 20, 'k', 'filled')
ylabel('Principal component')
title('Chlorophyta')
ylim([-13 15])

score_c = score; %Save this for later

% now do it all again for small phytoplankton group

rawvalues = table2array(small_phyto_Species(:,10:134));
n = height(rawvalues); 
offset = rawvalues+1; 
CLR = log(offset) - (1/n) .* sum(log(offset)); 

[wcoeff_picos,score,latent,tsquared,explained] = pca(CLR');

disp('Variance Explained:')
disp(explained(1:2))

subplot(4,2,1)
scatter(-score(:,1), score(:,2), 35, day(meta_noreps.sampledates, 'dayofyear'), 'filled')
xlabel('Principal component 1 (17.5%)')
ylabel('Principal component 2 (12.1%)')
title('Small Phytoplankton')
caxis([0 366])
ylim([-15 15])
xlim([-22 19])
box on 

subplot(4,1,2)
cla
plot(meta_noreps.sampledates(chrono_index), score(chrono_index,2), 'color', [.8 .8 .8], 'linewidth', 2)
hold on 
scatter(meta_noreps.sampledates(chrono_index), score(chrono_index,2), 20, [.8 .8 .8], 'filled')
plot(meta_noreps.sampledates(chrono_index), -score(chrono_index,1), 'color', [.2 .2 .2], 'linewidth', 2)
scatter(meta_noreps.sampledates(chrono_index), -score(chrono_index,1), 20, 'k', 'filled')
ylabel('Principal component')
title('Small Phytoplankton')
legend({'PC2'; ''; 'PC1'; ''})
ylim([-22 22])

% Finally, we want barplots with top contributors to each PC 
map = brewermap(12, 'paired');
temp = map(1,:);
map([1 4 5 8 9 12], :) = [];
map1 = brewermap(12, 'Set3');
map(4, :) = map1(end,:);
map(6, :) = temp;

% get top contributors to each PC for small phytos
[a, b] = sort(abs(wcoeff_picos(:,1)), 'descend'); 
subset = b(1:10); 
[a1, b1] = sort(-wcoeff_picos(subset), 'descend'); 

subplot(4,2,7)
cla
b = barh([1 2 3 4 5 6 7 8 9 10], -wcoeff_picos(subset(b1([1 2 3 4 5 6 7 8 9 10])))); 
hold on 
hold on
b.FaceColor = [.8 .8 .8]; 
order2 = [0 1 0 0 0 0 0 0 0 4]
for i = [2 10] 
b = barh(i, -wcoeff_picos(subset(b1(i)))); 
    hold on 
    b.FaceColor = map(order2(i),:); 
    b.FaceAlpha= .8
end

xlim([-.32 .32])
xlabel('Principal component 1 coefficient')
set(gca, 'Ydir', 'reverse')
yticks([1:10])
yticklabels(small_phyto_Species.Species(subset(b1)))
yticklabels({'\it Minutocellus polymorphus','\it Micromonas bravo B2','\it Chlorodendrales sp.','\it Gephyrocapsa oceanica','\it MOCH-3 sp.','\it Chrysophyceae Clade C', '\it Mantoniella squamata', '\it Pedinellales sp.','\it Dictyocha speculum','\it Phaeocystis pouchetii'});
set(gca,'TickLabelInterpreter','tex')
title('Small Phytoplankton')


% get top contributors to each PC for chlorophyta
[a, b] = sort(abs(wcoeff_chlorophytes(:,1)), 'descend'); 
subset = b(1:10); 
[a1, b1] = sort(wcoeff_chlorophytes(subset), 'descend'); 

subplot(4,2,8)
cla
b = barh([3 5 6 8 9 10], wcoeff_chlorophytes(subset(b1([3 5 6 8 9 10])))); 
hold on 
b.FaceColor = [.8 .8 .8]; 
order3 = [1 5 0 1 0 0 3];
for i = [1 2 4 7]
b = barh(i, wcoeff_chlorophytes(subset(b1(i)))); 
    hold on 
    b.FaceColor = map(order3(i),:); 
    b.FaceAlpha= .8;
end
xlim([-.32 .46])
xlabel('Principal component 1 coefficient')
set(gca, 'Ydir', 'reverse')
yticks([1:10])
yticklabels(chlorophyta_Species.Species(subset(b1)))
yticklabels({'\it Micromonas bravo B2','\it Ostreococcus Clade B','\it Chlorodendrales sp.','\it Micromonas pusilla','\it Pyramimonas australis','\it Pycnococcaceae sp.','\it Picochlorum sp.','\it Pterosperma sp.','\it Mantoniella squamata','\it Dolichomastigaceae-B sp.'});
set(gca,'TickLabelInterpreter','tex')
set(gca, 'Ydir', 'reverse')
title('Chlorophyta')

score_s = score; 

set(findall(gcf,'-property','FontSize'),'FontSize',15)

clear a1 a2 b1 g G h i ind2 ind2use subset explained order3 order1 orderr n offset order2 score rawvalues latent map1 a b tsquared 

%% Load in MVCO Flow Cytometry data for comparison 

% load('MVCO_timeseries.mat')
% 
% FCB_timeseries = timetable(); 
% 
% for yr = 2013:2021
%     path = [yr '\data\processed\grouped\'];
%     temp = load([path 'groupsum.mat']);
% 
%     ind2use = ~isnan(temp.cellresultsall(:,1));
% 
%     sampletime = datetime(temp.cellresultsall(ind2use,1), 'convertfrom', 'datenum');
%     euk_count = temp.cellNUMall(ind2use,4);
%     volume_analyzed = temp.cellresultsall(ind2use,3); 
% 
%     FCB_timeseries = [FCB_timeseries; timetable(sampletime, euk_count./volume_analyzed)];
% end
% FCB_timeseries = renamevars(FCB_timeseries,["Var1"],["euk_conc"]);
% 

load('FCB_timeseries.mat')

clear ind2use temp path euk_count volume_analyzed sampletime yr alleukconc allmatdate allsynconc 

%% Add daily average Euk concentration from FCB to metadata 
meta_noreps.FCM_euk_conc = nan(height(meta_noreps), 1); 

dailyAverage = retime(FCB_timeseries,'daily','mean');

for s = 1:height(meta_noreps)
    ind =find(dailyAverage.sampletime == meta_noreps.sampledates(s)); 
    if ~isempty(ind)
    meta_noreps.FCM_euk_conc(s) = dailyAverage.euk_conc(ind);
    end
end

figure
subplot(1,2,1)
scatter(meta_noreps.FCM_euk_conc, shanon_smallphyto, 50, table2array(sum(seqtab_noreps(:,19:143))), 'filled')
ylabel('Shannon diversity')
xlabel('Phytoeukaryote concentration (ml^{-1})')
title('Small Phytoplankton')
%set(gca,'ColorScale','log')
caxis([30000 300000]);

subplot(1,2,2)
scatter(meta_noreps.FCM_euk_conc, shanon_chlorophytes, 50, table2array(sum(seqtab_noreps(:,19:143))), 'filled')
xlabel('Phytoeukaryote concentration (ml^{-1})')
title("Chlorophytes")
%set(gca,'ColorScale','log')
caxis([30000 300000]);

colormap('viridis')
h = colorbar ;
h.Label.String = 'Sample depth (total reads)';

set(findall(gcf,'-property','FontSize'),'FontSize',20)

clear s ind temp_subset h latent map P R score explained tsquared ans 
%% Compare seasonality of concentration to seasonality of PCA

%First plot top row of PCA scores over time 
map2 = viridis(12);
figure
subplot(2,2,2)
[~, chrono_index] = sort(meta_noreps.sampledates); %get chronological index so we can plot connected dots 
for yr = 2013:2021
 yrind = find(year(meta_noreps.sampledates(chrono_index)) == yr);
 hold on 
    scatter(day(meta_noreps.sampledates(chrono_index(yrind)), 'dayofyear'), score_c(chrono_index(yrind),1), 20, map2(yr-2012, :), 'filled')
    plot(day(meta_noreps.sampledates(chrono_index(yrind)), 'dayofyear'), score_c(chrono_index(yrind),1), 'linewidth', 2, 'color', map2(yr-2012, :))
end
xlim([0 366])

subplot(2,2,1)
for yr = 2013:2021
 yrind = find(year(meta_noreps.sampledates(chrono_index)) == yr);
 hold on 
    scatter(meta_noreps.sampledates(chrono_index(yrind)), score_c(chrono_index(yrind),1), 20, map2(yr-2012, :), 'filled')
    plot(meta_noreps.sampledates(chrono_index(yrind)), score_c(chrono_index(yrind),1), 'linewidth', 2, 'color', map2(yr-2012, :))
end

ylabel({'Chlorophyta'; 'Principle Component 1'}, 'FontSize', 20)



%Now do second row of concentrations over time. 

% Load in Discrete data to fill gaps in time serise

%discretes = readtable('\\MVCO_FCMdiscrete\nes-lter-fcm-discrete-mvco.csv');
%discretes = discretes(discretes.depth_m < 10, :); %only use surface samples
%sampletime = discretes.date_sampled;
%euk_conc = discretes.redeuk_leq_20um_cells_per_ml;
%FCM_timeseries = timetable(sampletime, euk_conc);
% save('FCM_timeseries.mat', 'FCM_Timeseries')


load('FCM_timeseries.mat')

subplot(2,2,4)

for yr = 2013:2021
 yrind = find(year(dailyAverage.sampletime) == yr);
 hold on 
   scatter(day(dailyAverage.sampletime(yrind), 'dayofyear'), dailyAverage.euk_conc(yrind), 20, map2(yr-2012, :), 'filled')
    plot(day(dailyAverage.sampletime(yrind), 'dayofyear'), dailyAverage.euk_conc(yrind), 'linewidth', 2, 'color', map2(yr-2012, :))
end

subplot(2,2,3)
for yr = 2013:2021
 yrind = find(year(dailyAverage.sampletime) == yr) ;
 hold on 
    scatter(dailyAverage.sampletime(yrind), dailyAverage.euk_conc(yrind), 20, map2(yr-2012, :), 'filled')
    plot(dailyAverage.sampletime(yrind), dailyAverage.euk_conc(yrind), 'linewidth', 2, 'color', map2(yr-2012, :))
    
 yrind = find(year(FCM_timeseries.sampletime) == yr) ;
 scatter(FCM_timeseries.sampletime(yrind), FCM_timeseries.euk_conc(yrind), 20, map2(yr-2012, :), '^')

end

ylabel({'Phytoeukaryote concentration'; 'measured by FlowCytobot (cells ml^{-1})'},'FontSize',20)

legend({'Daily average', '', 'Discrete samples'})
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    ylim([0 10.9e4])

subplot(2,2,4)
xlim([0 366])
xlabel('Day of year', 'FontSize', 20)

clear yrind ind j l 


%% Bring in environmental data for supplementary figure 
% 
%load('\\EnvironmentalData\MVCO_Environmental_Tables.mat')
%nutrients = readtable('\\EnvironmentalData\nes-lter-nutrient-mvco.csv');
%S1 = readtable('2020_Radiation_asit.txt');
%S2 = readtable('2021_Radiation_asit.txt');
%save('Environmental_Data_FigureS3.mat', 'MVCO_Env_Table', 'nutrients', 'S1', 'S2', 'MVCO_Daily')

load('Environmental_Data_FigureS3.mat'); 

[groupnum, id] = findgroups(floor(datenum(MVCO_Env_Table.time_local)));
AvgSalinity = splitapply(@nanmean, MVCO_Env_Table.salinity_beam, groupnum);
[~,iA] = intersect(datenum(MVCO_Daily.days),id);
MVCO_Daily.AvgSalinity(iA) = AvgSalinity;

PC1 = -score_s(:,1);

figure
subplot(5,1,1)
% want a climatology behind each plot in grey
meta_noreps.doy = day(meta_noreps.sampledates, 'dayofyear'); 
meta_noreps.pc1 = PC1; 
climat = groupsummary(meta_noreps, ["doy"], "mean", "pc1");
hold on 
for y = 2013:2021; 
    timevec = datenum(y, 1, 1) + climat.doy ;
    plot(datetime(timevec, 'convertfrom', 'datenum'), climat.mean_pc1, 'color', [.8 .8 .8], 'linewidth', 1.2)
end
scatter(meta_noreps.sampledates, PC1, 10, 'k', 'filled')
ylabel('PC1')
a = gca; 
xlimits = a.XLim;


subplot(5,1,2)
% want a climatology behind each plot in grey
MVCO_Daily.doy = day(MVCO_Daily.days, 'dayofyear'); 
climat = groupsummary(MVCO_Daily, ["doy"], "mean", "Beam_temperature_corrected");
hold on 
for y = 2013:2021; 
    timevec = datenum(y, 1, 1) + climat.doy ;
    plot(datetime(timevec, 'convertfrom', 'datenum'), climat.mean_Beam_temperature_corrected, 'color', [.8 .8 .8], 'linewidth', 1.2)
end
scatter(MVCO_Daily.days, MVCO_Daily.Beam_temperature_corrected, 10, 'k', 'filled')
ylabel(['Daily Average Temperature (' char(176) 'C)'])
xlim(xlimits)

MVCO_Env_Table.doy = day(MVCO_Env_Table.time_local, 'dayofyear'); 

subplot(5,1,3)
climat = groupsummary(MVCO_Env_Table(year(MVCO_Env_Table.time_local)>2012, :), ["doy"], "mean", "salinity_beam");
hold on 
for y = 2013:2021; 
    timevec = datenum(y, 1, 1) + climat.doy ;
    plot(datetime(timevec, 'convertfrom', 'datenum'), climat.mean_salinity_beam, 'color', [.8 .8 .8], 'linewidth', 1.2)
end
scatter(MVCO_Daily.days, MVCO_Daily.AvgSalinity, 10, 'k', 'filled')
ylabel("Daily Average Salinity (PSU)")
ylim([29.5 32.5])
xlim(xlimits)

subplot(5,1,4)
nutrients = nutrients(nutrients.depth<10,:); 
nutrients = nutrients(year(nutrients.date_time_utc) > 2012,:); 
nutrients.doy = day(nutrients.date_time_utc, 'dayofyear'); 
climat = groupsummary(nutrients, ["doy"], "mean", "phosphate");
hold on 
for y = 2013:2021; 
    timevec = datenum(y, 1, 1) + climat.doy ;
    plot(datetime(timevec, 'convertfrom', 'datenum'), climat.mean_phosphate, 'color', [.8 .8 .8], 'linewidth', 1.2)
end
scatter(nutrients.date_time_utc, nutrients.phosphate, 10, 'k', 'filled')
xlim(xlimits)
ylabel("Phosphate (/mu M)")


subplot(5,1,5)
climat = groupsummary(MVCO_Daily, ["doy"], "mean", "AvgSolar");
hold on 
for y = 2013:2021; 
    timevec = datenum(y, 1, 1) + climat.doy ;
    plot(datetime(timevec, 'convertfrom', 'datenum'), climat.mean_AvgSolar, 'color', [.8 .8 .8], 'linewidth', 1.2)
end
xlim(xlimits)

MVCO_Env_Table.month = month(MVCO_Env_Table.time_local);
MVCO_Env_Table.year = year(MVCO_Env_Table.time_local);
means = groupsummary(MVCO_Env_Table, ["month" "year"], "mean", "solar");
scatter(datetime(means.year, means.month, 1), means.mean_solar, 15, 'k', 'filled')
ylabel(['Monthly Average Solar Radiation (W m^{-2})'])
xlim(xlimits)


S = [S1; S2]; 
S.sampledate = datetime(S.Var1(:), 'InputFormat', 'yyyy-mm-ddTHH:MM:ss');
S.month = month(S.sampledate); 
S.year = year(S.sampledate);
means = groupsummary(S, ["month" "year"], "mean", "Var2");
hold on 
scatter(datetime(means.year, means.month, 1), means.mean_Var2, 15, 'k', 'filled')
xlim(xlimits)

clear climat a map2 means iA id S S1 S2 timevec y xlimits AvgSalinity groupnum


%% Make Figure 3 with 5 panels. 
% Load division rate data from Stevens et al. 
load('Data_Fig4.mat')

figure

subplot(5,1,2)
scatter(datetime(day(meta_noreps.sampledates, 'dayofyear'), 'convertfrom', 'datenum'), -score_s(:,1), 20, [.3 .3 .3], 'filled', 'MarkerFaceAlpha', .8)

PC1_climat = nan(1,12);
PC1_std = nan(1,12); 
for i = 1:12
    ind = find(month(meta_noreps.sampledates) == i);
    PC1_climat(i) = nanmean(-score_s(ind,1)); 
    PC1_std(i) = nanstd(-score_s(ind,1)); 
end

 
hold on 
h = area(datetime(0000, 1:12, 15), [PC1_climat-PC1_std; 2.*PC1_std]'); 
h(1).FaceAlpha = 0;
h(2).FaceAlpha = .15; 
h(2).FaceColor = 'k'; 
h(2).EdgeColor = 'none'; 
h(1).EdgeColor = 'none'; 
hold on 
h(1).BaseLine.LineStyle='none'

plot(datetime(0000, 1:12, 15), PC1_climat, 'k', 'linewidth', 2)

ylabel({'Small Phytoplankton'; 'principle component 1'})
a = gca; 
xlimits = a.XLim; 



subplot(5,1,4)
scatter(datetime(day(meta_noreps.sampledates, 'dayofyear'), 'convertfrom', 'datenum'), score_c(:,1), 20, [.3 .3 .3], 'filled', 'MarkerFaceAlpha', .8)

PC1_climat = nan(1,12);
PC1_std = nan(1,12); 
for i = 1:12
    ind = find(month(meta_noreps.sampledates) == i);
    PC1_climat(i) = nanmean(score_c(ind,1)); 
    PC1_std(i) = nanstd(score_c(ind,1)); 
end

 
hold on 
h = area(datetime(0000, 1:12, 15), [PC1_climat-PC1_std; 2.*PC1_std]'); 
h(1).FaceAlpha = 0;
h(2).FaceAlpha = .15; 
h(2).FaceColor = 'k'; 
h(2).EdgeColor = 'none'; 
h(1).EdgeColor = 'none'; 
hold on 
h(1).BaseLine.LineStyle='none'

plot(datetime(0000, 1:12, 15), PC1_climat, 'k', 'linewidth', 2)

ylabel({'Chlorophyta'; 'principle component 1'})
a = gca; 
xlimits = a.XLim; 

subplot(5,1,3)
 
FCB_timeseries = FCB_timeseries(year(FCB_timeseries.sampletime) > 2012,:);
FCB_climat = nan(1,12); 
FCB_std = nan(1,12)
for i = 1:12
    ind = find(month(FCB_timeseries.sampletime) == i);
    FCB_climat(i) = nanmean(FCB_timeseries.euk_conc(ind)); 
    FCB_std(i) = nanstd(FCB_timeseries.euk_conc(ind)); 
end
hold on 
h = area(datetime(0000, 1:12, 15), [FCB_climat-FCB_std; 2.*FCB_std]'); 
h(1).FaceAlpha = 0;
h(2).FaceAlpha = .15; 
h(2).FaceColor = 'k'; 
h(2).EdgeColor = 'none'; 
h(1).EdgeColor = 'none'; 
hold on 
h(1).BaseLine.LineStyle='none'

dailyAverage = dailyAverage(year(dailyAverage.sampletime) > 2012,:);

plot(datetime(0000, 1:12, 15), FCB_climat, 'k', 'linewidth', 2)
scatter(day(dailyAverage.sampletime, 'dayofyear'), dailyAverage.euk_conc,  20, [.3 .3 .3], 'filled', 'MarkerFaceAlpha', .8)
ylabel({'Phytoeukaryotes'; '(cells ml^{-1})'},'FontSize',20)
xlim(xlimits)

subplot(5,1,1)
scatter(datetime(1:366, 'convertfrom', 'datenum'), nanmean(Daily_temp'),  20, [.3 .3 .3], 'filled', 'MarkerFaceAlpha', .8)

temp_climat = nan(1,12);
temp_std = nan(1,12); 
monthvec = month(datetime(1:366, 'convertfrom', 'datenum')); 

for i = 1:12
    ind = find(monthvec == i);
    allmonthstemp = Daily_temp(ind,:);
    temp_climat(i) = nanmean(allmonthstemp(:)); 
    temp_std(i) = nanstd(allmonthstemp(:)); 
end

hold on 
h = area(datetime(0000, 1:12, 15), [temp_climat-temp_std; 2.*temp_std]'); 
h(1).FaceAlpha = 0;
h(2).FaceAlpha = .15; 
h(2).FaceColor = 'k'; 
h(2).EdgeColor = 'none'; 
h(1).EdgeColor = 'none'; 
hold on 
h(1).BaseLine.LineStyle='none'

plot(datetime(0000, 1:12, 15), temp_climat, 'k', 'linewidth', 2)

ylabel(['Temperature (' char(176) 'C)'])

    set(findall(gcf,'-property','FontSize'),'FontSize',16)

subplot(5,1,5)
scatter(datetime(1:366, 'convertfrom', 'datenum'), nanmean(daily_euk_divrate'),  20, [.3 .3 .3], 'filled', 'MarkerFaceAlpha', .8)

dr_climat = nan(1,12);
dr_std = nan(1,12); 
monthvec = month(datetime(1:366, 'convertfrom', 'datenum')); 

for i = 1:12
    ind = find(monthvec == i);
    allmonthsdivrates = daily_euk_divrate(ind,:);
    dr_climat(i) = nanmean(allmonthsdivrates(:)); 
    dr_std(i) = nanstd(allmonthsdivrates(:)); 
end

hold on 
h = area(datetime(0000, 1:12, 15), [dr_climat-dr_std; 2.*dr_std]'); 
h(1).FaceAlpha = 0;
h(2).FaceAlpha = .15; 
h(2).FaceColor = 'k'; 
h(2).EdgeColor = 'none'; 
h(1).EdgeColor = 'none'; 
hold on 
h(1).BaseLine.LineStyle='none'

plot(datetime(0000, 1:12, 15), dr_climat, 'k', 'linewidth', 2)

ylabel({'Phytoeukaryote'; 'daily division rate (d^{-1})'})
ylim([0 3])

    set(findall(gcf,'-property','FontSize'),'FontSize',16)

    clear PC1_climat PC1_std R p P dr_climat FCB_climat FCB_std dr_std s score temp_climat temp_std temp_subset allmonthsdivrates allmonthstemp monthvec mse_temp mse_doy 


%%  Supplementary Figure S4
% want a supplementary figure with correlations between PC1 and temp, PC 1 and conc, 

figure 
subplot(3,2,1)
scatter(meta_noreps.temperature, -score_s(:,1))
[R, P] = corrcoef(meta_noreps.temperature, -score_s(:,1), 'rows', 'complete')
text(1, 15, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
xlim([0 25])
ylabel('Small Phytoplankton PC1')
xlabel(['Temperature (' char(176) 'C)'])
subplot(3,2,2)
scatter(meta_noreps.temperature, score_c(:,1))
[R, P] = corrcoef(meta_noreps.temperature, score_c(:,1), 'rows', 'complete');
text(1, 10, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
xlim([0 25])
xlabel(['Temperature (' char(176) 'C)'])
ylabel('Chlorophyta PC1')

subplot(3,2,3)
scatter(meta_noreps.FCM_euk_conc, -score_s(:,1))
[R, P] = corrcoef(meta_noreps.FCM_euk_conc, -score_s(:,1), 'rows', 'complete')
ylim([-20 22])
text(4200, 17, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
ylabel('Small Phytoplankton PC1')
xlabel('Phytoeukaryote Concentration (ml^{-1})')

subplot(3,2,4)
scatter(meta_noreps.FCM_euk_conc, score_c(:,1))
[R, P] = corrcoef(meta_noreps.FCM_euk_conc, score_c(:,1), 'rows', 'complete')
text(4200, 15, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
ylim([-10 18])
xlabel('Phytoeukaryote Concentration (ml^{-1})')
ylabel('Chlorophyta PC1')

subplot(3,2,5)
scatter(meta_noreps.phosphate, -score_s(:,1))
[R, P] = corrcoef(meta_noreps.phosphate, -score_s(:,1), 'rows', 'complete')
ylim([-20 22])
text(.6, 17, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
ylabel('Small Phytoplankton PC1')
xlabel('Phosphate (\mu M)')

subplot(3,2,6)
scatter(meta_noreps.phosphate, score_c(:,1))
[R, P] = corrcoef(meta_noreps.phosphate, score_c(:,1), 'rows', 'complete')
text(.6, 15, {['R = ' num2str(R(1,2))]; ['p = ' num2str(P(1,2))]})
ylim([-10 18])
xlabel('Phosphate (\mu M)')
ylabel('Chlorophyta PC1')


%%
% Compare predictive power of day of year vs temperature

% Polynomial of degree 3 fit to dayof year has adjusted R squared 0.53 and
% RMSE 3.58

%polynomial of degree 3 fit to temperature data has adjusted R squared of
%0.566 and RMS 3.56

%% Formal test of null hypothesis that temperature and day of year are equally good predictors

X_doy = day(meta_noreps.sampledates, 'dayofyear'); 
X_temp = meta_noreps.temperature; 
y = score_s(:,1);

X_doy(isnan(X_temp)) = []; 
y(isnan(X_temp)) = []; 
X_temp(isnan(X_temp)) = []; 

degree = 3; % best for nearly sinusoidal functions

K = 10;
cv = cvpartition(length(y), 'KFold', K);

mse_temp = zeros(K,1);
mse_doy = zeros(K,1);

for i = 1:K
    trainIdx = training(cv, i);
    testIdx = test(cv, i);

    % Training data
    X_temp_train = X_temp(trainIdx);
    X_doy_train = X_doy(trainIdx);
    y_train = y(trainIdx);

    % Test data
    X_temp_test = X_temp(testIdx);
    X_doy_test = X_doy(testIdx);
    y_test = y(testIdx);

    % Fit polynomial models
    p_temp = polyfit(X_temp_train, y_train, degree);
    p_doy = polyfit(X_doy_train, y_train, degree);

    % Predict
    y_pred_temp = polyval(p_temp, X_temp_test);
    y_pred_doy = polyval(p_doy, X_doy_test);

    % Compute mean squared error
    mse_temp(i) = mean((y_test - y_pred_temp).^2);
    mse_doy(i) = mean((y_test - y_pred_doy).^2);

    subplot(1,2,1)
    hold on 
    scatter(y_pred_temp, y_test)
    
    subplot(1,2,2)
    hold on 
    scatter(y_pred_doy, y_test)

end

% Paired t-test
[h, p, ci, stats] = ttest(mse_temp, mse_doy);

% Display results
fprintf('Mean MSE (temperature): %.4f\n', mean(mse_temp));
fprintf('Mean MSE (day-of-year): %.4f\n', mean(mse_doy));
fprintf('Paired t-test p-value: %.4f\n', p);

clear i p P p_day p_temp R PC1 ci cv a h ind degree p_doy K y_pred_temp y_test X_doy X_doy_test X_doy_train stats y testIdx trainIdx mse_temp mse_doy  xlimits X_temp_test X_temp X_temp_train y_pred_doy y_train yr
%% Generate Phaeosystis supplement figure 

temp = small_phyto_noreps(small_phyto_noreps.Species_boot>80, :);
small_phyto_Species = grpstats(temp, {'Kingdom'; 'Supergroup'; 'Division'; 'Class'; 'Order'; 'Family'; 'Genus'; 'Species'}, {'sum'}, 'DataVars', temp.Properties.VariableNames(19:143));
small_phyto_Species.total = sum(small_phyto_Species{:,10:134},2);
small_phyto_Species = sortrows(small_phyto_Species, 'total', 'descend');

%have to reorder chronologically
[~, chrono_index] = sort(meta_noreps.sampledates); %get chronological index so we can plot connected dots 

%Now use these tables to calculate CLR 
rawvalues = table2array(small_phyto_Species(:,10:134));
n = height(rawvalues); 
offset = rawvalues+1; 
CLR = log(offset) - (1/n) .* sum(log(offset)); 
CLR_chrono = CLR(:,chrono_index);

phaeo_ind = find(contains(small_phyto_Species.Species, "Phaeocystis"));
phaeo_sum = sum(CLR_chrono(phaeo_ind, :));

figure
yyaxis left
scatter(meta_noreps.sampledates(chrono_index), phaeo_sum, 20, 'b', 'filled')
hold on 
plot(meta_noreps.sampledates(chrono_index), phaeo_sum, 'b')
datetick
b = gca; 
b.YColor = 'b';
ylim([-8 16])
ylabel('Sum CLR Phaeocystis')

%Load FCB data
%load('\\FCB_compiledC_tables.mat')
%save('TTsumnum.mat', 'TTsumnum')
load('TTsumnum.mat')

yyaxis right

bigeukconc = TTsumnum.Eukleq10microns./TTsumnum.ml - TTsumnum.Eukleq3microns./TTsumnum.ml; 
[G, i]  = findgroups(dateshift(TTsumnum.Time, 'start', 'day'))
M = splitapply(@nanmean, bigeukconc, G)
scatter(i, M, 10, 'k', 'filled')
hold on 
%plot(i, M, 'k')
xlim([min(meta_noreps.sampledates), max(meta_noreps.sampledates)])
ylim([-1000 6000])
ylabel('FCB eukaryotes 3-10 (\mum^{3})')

a = gca; 
a.YColor = 'k';

set(findall(gcf,'-property','FontSize'),'FontSize',20)


%% function just cuts table to chlorophtyes, kinda silly

function cutable = cut_to_chlorophytes(speciestable)

cutable = speciestable; 
cutable(~strcmp(cutable.Division, 'Chlorophyta'), :) = [];

end