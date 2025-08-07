% Use this script to visualize results of topic selection from ROST CLI

%% Load in data and sort into necessary parts

%set location where outputs of rostcli are for chlorophyta 
foldername = 'chlorophyte_topicmodel' ; 

%add path to brewermap so we can have pretty colors
addpath('\\DrosteEffect-BrewerMap-3.2.5.0')


% Load topic model results for chlorophyta 
topic_model_c = readtable(strcat( foldername, "/topicmodel.csv"));
topic_model_c = topic_model_c{:,:};
topic_maxl_c = readmatrix(strcat(foldername, "/topics.maxlikelihood.csv"));
perp_c = readmatrix(strcat(foldername, "/perplexity.iter.csv"));

Res_c = topic_maxl_c(:,2:end); 
timestamps = topic_maxl_c(:,1); 


N_topics = 2; %This is the value of the paramter k from topics.refine command 
% Make a nice table with # of observations from each topic within each
% sample
Hist_c = []; 
for i = 1:size(Res_c, 1)   
    X = histcounts(Res_c(i,:), 0:N_topics); 
    Hist_c = [Hist_c; X]; 
    clear X
end


%set location where outputs of rostcli are for small phytoplanton 
foldername = 'small_phyto_topicmodel' ; 


% Load topic model results of small phytoplankton 
topic_model_p = readtable(strcat( foldername, "/topicmodel.csv"));
topic_model_p = topic_model_p{:,:};
topic_maxl_p = readmatrix(strcat(foldername, "/topics.maxlikelihood.csv"));
perp_p = readmatrix(strcat(foldername, "/perplexity.iter.csv"));

Res_p = topic_maxl_p(:,2:end); 

timestamps_p = topic_maxl_p(:,1); 

N_topics = 2; %This is the value of the paramter k from topics.refine command 
% Make a nice table with # of observations from each topic within each
% sample
Hist_p = []; 
for i = 1:size(Res_p, 1)   
    X = histcounts(Res_p(i,:), 0:N_topics); 
    Hist_p = [Hist_p; X]; 
    clear X
end


% We also want to load metadata  & sort it chronologically as we did before applying ROST

% Complete metadata in original order
   sorted_meta = readtable("MVCO_metadata_AllSequenceSamples.csv");
    
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
%now sort chronologically 
[a, b]= sort(meta_noreps.sampledates); 
chrono_meta_noreps = meta_noreps(b,:);

%% Make stacked barplot of topic hists 

figure
subplot(2,1,2)

Scaled_hist = Hist_c'./sum(Hist_c'); 
dates = chrono_meta_noreps.sampledates; 


map = brewermap(11,'spectral');
map(3:6,:) = []; 
map = [map(3:end, :); map(1:2, :)];
map = map(end:-1:1, :);


Scaled_hist = [Scaled_hist(2,:); Scaled_hist(1,:)] %flip them because we want spring community to be blue and on top

a = bar(timestamps+min(dates), Scaled_hist', 20, 'stacked'); 
for i = 1:N_topics;
    a(i).FaceColor = map(i,:); 
    a(i).FaceAlpha= (0.75); 
    a(i).EdgeColor = 'none'; 
end
ylabel('Proportion of reads')
ylim([0 1])

a(1).FaceColor = map(2,:)
a(2).FaceColor = map(5,:)
title('Chlorophyta')


subplot(2,1,1)

Scaled_hist = Hist_p'./sum(Hist_p'); 
Scaled_hist = [Scaled_hist(2,:); Scaled_hist(1,:)]

a = bar(timestamps_p+min(dates), Scaled_hist', 20, 'stacked'); 
for i = 1:N_topics;
    a(i).FaceColor = map(i,:); 
    a(i).FaceAlpha= (0.75); 
    a(i).EdgeColor = 'none'; 
end
ylabel('Proportion of reads')
ylim([0 1])
legend
a(1).FaceColor = map(2,:)
a(2).FaceColor = map(5,:)
title('Small Phytoplankton')
legend({'Co-occurrence Community 2'; 'Co-occurrence Community 1'})


%% Now look at taxa within these communities 

figure

map = brewermap(11,'spectral');
map(3:6,:) = []; 
map = [map(3:end, :); map(1:2, :)];
map = map(end:-1:1, :);


subplot(2,1,1)

Scaled_hist = Hist_p'./sum(Hist_p'); 
Scaled_hist = [Scaled_hist(2,:); Scaled_hist(1,:)]

a = bar(timestamps_p+min(dates), Scaled_hist', 20, 'stacked'); 
for i = 1:N_topics;
    a(i).FaceColor = map(i,:); 
    a(i).FaceAlpha= (0.75); 
    a(i).EdgeColor = 'none'; 
end
ylabel('Proportion of reads')
ylim([0 1])
legend
a(1).FaceColor = map(2,:)
a(2).FaceColor = map(5,:)
title('Small Phytoplankton')
legend({'Co-occurrence Community 2'; 'Co-occurrence Community 1'})


%load word dictionary from part 4-1
load('small_phyto_topicmodel\SmallPhyto_wordlist.mat')

[~, maxorder] = sort(topic_model_p(2,:), 'descend'); 

scaled_topic_model = topic_model_p./sum(topic_model_p,2);

% Color stuff
map = brewermap(12, 'paired');
temp = map(1,:);
map([1 4 5 8 9 12], :) = [];
map1 = brewermap(12, 'Set3');
map(4, :) = map1(end,:);
map(6, :) = temp;

subplot(2,2,4)
b = bar([1:10], scaled_topic_model(2, maxorder([1:10])))
hold on
%b.BarWidth = .12
b.FaceColor = [.8 .8 .8]; 
order = [3 0 0 0 2 1 1 0 5 5 0];
for i = [1 5 6 7 9 10 ]%[1 2 4 5 6 7 8 9] 
    b = bar(i, scaled_topic_model(2, maxorder(i)))
    hold on 
    b.FaceColor = map(order(i),:); 
    b.FaceAlpha= .8
end

% 1 is micormoas, 2 is bathyococus, 3 is picochlorum, 6
% ostreococcus, 5 phaeo

xticks([1:10])
labels = string(word_dict(maxorder(1:10))); 
%labels = {'\it Picochlorum sp.','\it Chlorodendrales sp.','\it MOCH-3 sp.', '\it Minutocellus polymorphus', '\it Bathycoccus prasinos', '\it Micromonas bravo B2', '\it Micromonas commoda A2', '\it Gephyrocapsa oceanica', '\it Ostreococcus lucimarinus', '\it Ostreococcus Clade B'};

xticklabels(labels)
title('Co-occurrence Community 2')
ax = gca; 
ax.TickLabelInterpreter = 'tex';
ylabel('Probability')


[~, maxorder] = sort(topic_model_p(1,:), 'descend'); 

scaled_topic_model = topic_model_p./sum(topic_model_p,2);

subplot(2,2,3)
b = bar([1:10], scaled_topic_model(1,maxorder([1:10])))
%b.FaceColor = map(2, :);
b.FaceColor = [.8 .8 .8]; 
order = [2 1 0 4 0 0 0 0 0 1]% [2 1 4 9 1]
for i = [1 2 4 10]
    hold on 
    b = bar(i, scaled_topic_model(1, maxorder(i)))
    b.FaceColor = map(order(i),:); 
    b.FaceAlpha= .8
end
ylabel('Probability')
labels = string(word_dict(maxorder(1:10))); 
%labels = {'\it Bathycoccus prasinos', '\it Micromonas commoda A2', '\it Aureococcus anophagefferens','\it Phaeocystis pouchetii','\it Plagioselmis prolonga', '\it Pterosperma sp.', '\it Falcomonas daucoides',  '\it Chrysophyceae Clade-C', '\it Pedinellales sp.', '\it Micromonas pusilla'};

xticks([1:10])
xticklabels(labels)
title('Co-occurrence Community 1')
ax = gca; 
ax.TickLabelInterpreter = 'Tex';

set(findall(gcf,'-property','FontSize'),'FontSize',16)

%% 

% and now for chlorophyta % 

load('chlorophyte_topicmodel\Chlorophyta_wordlist.mat')

[~, maxorder] = sort(topic_model_c(1,:), 'descend'); 
scaled_topic_model = topic_model_c./sum(topic_model_c,2);

figure

map = brewermap(11,'spectral');
map(3:6,:) = []; 
map = [map(3:end, :); map(1:2, :)];
map = map(end:-1:1, :);


Scaled_hist = Hist_c'./sum(Hist_c'); 
Scaled_hist = [Scaled_hist(2,:); Scaled_hist(1,:)]


subplot(2,1,1)
a = bar(timestamps+min(dates), Scaled_hist', 20, 'stacked'); 
for i = 1:N_topics;
    a(i).FaceColor = map(i,:); 
    a(i).FaceAlpha= (0.75); 
    a(i).EdgeColor = 'none'; 
end
ylabel('Proportion of reads','FontSize',20);
ylim([0 1])

a(1).FaceColor = map(2,:)
a(2).FaceColor = map(5,:)

legend1 = legend({'Co-occurrence Community 2'; 'Co-occurrence Community 1'});
set(legend1,...
    'Position',[0.727458468931138 0.451203424788978 0.145167229392887 0.103272931357106],...
    'FontSize',16);


title('Chlorophyta')

% Color stuff
map = brewermap(12, 'paired');
temp = map(1,:);
map([1 4 5 8 9 12], :) = [];
map1 = brewermap(12, 'Set3');
map(4, :) = map1(end,:);
map(6, :) = temp;

subplot(2,2,3)
b = bar([1:10], scaled_topic_model(1,maxorder([1:10])))
%b.FaceColor = map(2, :); 
hold on 
b.FaceColor = [.8 .8 .8]; 
order = [2 1 0 1 5 0 1 0 0 5]%[2 1 9 1 5 1 5];
for i = [1 2 4 5 7 10]
    b = bar(i, scaled_topic_model(1, maxorder(i)))
    b.FaceColor = map(order(i),:); 
    b.FaceAlpha= .8
end
xticks([1:10])
labels = string(word_dict(maxorder(1:10))); 
%labels4 = {'\it Bathycoccus prasinos','\it Micromonas commoda A2','\it Pterosperma sp.','\it Micromonas pusilla','\it Ostreococcus lucimarinus','\it Mamiella gilva','\it Micromonas bravo B2','\it Dolichomastigaceae B sp.','\it Pyramimonas australis','\it Ostreococcus Clade B'};


xticklabels(labels4)
title('Co-occurrence Community 1')
ax = gca; 
ax.TickLabelInterpreter = 'tex';

% 9 is pterosperma, 1 is micormoas, 2 is bathyococus, 3 is picochlorum, 5 ostreococcus, 6 is chlorodendrales 
% 4 is phaeo

[~, maxorder] = sort(topic_model_c(2,:), 'descend'); 

scaled_topic_model = topic_model_c./sum(topic_model_c,2);

subplot(2,2,4)
b = bar([1:10], scaled_topic_model(2,maxorder([1:10])))
hold  on 
b.FaceColor = [.8 .8 .8]; 
order = [3 0 1 0 5 0 0 0 1 0]%[3 6 1 5 1];
for i = [1 3 5 9]
    b = bar(i, scaled_topic_model(2, maxorder(i)))
    hold on 
    b.FaceColor = map(order(i),:); 
    b.FaceAlpha= .8
end
xticks([1:10])
labels = string(word_dict(maxorder(1:10))); 

%labels3 = {'\it Picochlorum sp.','\it Chlorodendrales sp.','\it Micromonas bravo B2','\it Pycnococcaceae sp.','\it Ostreococcus Clade B','\it Chlamydomonas parkeae','\it Pyramimonas australis','\it Chloroidium ellipsoideum','\it Micromonas pusilla','\it Trebouxia sp.'};

xticklabels(labels3)

title('Co-occurrence Community 2','FontSize',16);

ax = gca; 
ax.TickLabelInterpreter = 'tex';


set(findall(gcf,'-property','FontSize'),'FontSize',16)

%make y labels bigger
subplot(2,1,1)
ylabel('Proportion of reads','FontSize',20);
subplot(2,2,3)
ylabel('Probability','FontSize',20);
