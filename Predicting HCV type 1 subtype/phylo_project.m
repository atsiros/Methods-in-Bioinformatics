%% Load in sequences data file
filename = fopen('/Users/atsiros/full_seqs.csv');
cols = textscan(filename,'%f %s %s %s','delimiter',',');
fclose(filename);

%% Define column
sequences = cols{3};

%% Split into appropriate years
seq99 = sequences(42:91);
seq00 = sequences(92:107);
seq01 = sequences(108:112);
seq03 = sequences(113:143);
seq06 = sequences(147:153);
seq07 = sequences(154:310);
seq08 = sequences(311:349);
seq10 = sequences(350:358);

%% Create multiple sequence alignment and define length
MSA99 = multialign(seq99);
MSA00 = multialign(seq00);
MSA01 = multialign(seq01);
MSA03 = multialign(seq03);
MSA06 = multialign(seq06);
MSA07 = multialign(seq07);
MSA08 = multialign(seq08);
MSA10 = multialign(seq10);

%% Build PHMM model
model99 = hmmprofstruct(length(MSA99));
model99 = hmmprofestimate(model99,MSA99);

model00 = hmmprofstruct(length(MSA00));
model00 = hmmprofestimate(model00,MSA00);

model01 = hmmprofstruct(length(MSA01));
model01 = hmmprofestimate(model01,MSA01);

model03 = hmmprofstruct(length(MSA03));
model03 = hmmprofestimate(model03,MSA03);

model06 = hmmprofstruct(length(MSA06));
model06 = hmmprofestimate(model06,MSA06);

model07 = hmmprofstruct(length(MSA07));
model07 = hmmprofestimate(model07,MSA07);

model08 = hmmprofstruct(length(MSA08));
model08 = hmmprofestimate(model08,MSA08);

model10 = hmmprofstruct(length(MSA10));
model10 = hmmprofestimate(model10,MSA10);

%% Calculating scores from the models
scores99 = [];

for x = 1:length(seq99)
    score = hmmprofalign(model99,seq99{x});
    scores99 = [scores99; score];
end

scores00 = [];

for x = 1:length(seq00)
    score = hmmprofalign(model00,seq00{x});
    scores00 = [scores00; score];
end

scores01 = [];

for x = 1:length(seq01)
    score = hmmprofalign(model01,seq01{x});
    scores01 = [scores01; score];
end

scores03 = [];

for x = 1:length(seq03)
    score = hmmprofalign(model03,seq03{x});
    scores03 = [scores03; score];
end

scores06 = [];

for x = 1:length(seq06)
    score = hmmprofalign(model06,seq06{x});
    scores06 = [scores06; score];
end

scores07 = [];

for x = 1:length(seq07)
    score = hmmprofalign(model07,seq07{x});
    scores07 = [scores07; score];
end

scores08 = [];

for x = 1:length(seq08)
    score = hmmprofalign(model08,seq08{x});
    scores08 = [scores08; score];
end

scores10 = [];

for x = 1:length(seq10)
    score = hmmprofalign(model10,seq10{x});
    scores10 = [scores10; score];
end

%%  Tree preparation
[max99,i99] = max(scores99);
[max00,i00] = max(scores00);
[max01,i01] = max(scores01);
[max03,i03] = max(scores03);
[max06,i06] = max(scores06);
[max07,i07] = max(scores07);
[max08,i08] = max(scores08);
[max10,i10] = max(scores10);

seqHi = [seq99(i99) seq00(i00) seq01(i01) seq03(i03) seq06(i06) seq07(i07) seq08(i08) seq10(i10)];

distances = seqpdist(seqHi,'Method','Jukes-Cantor');

names = ['99';'00';'01';'03';'06';'07';'08';'10'];
names = cellstr(names);

tree = seqlinkage(distances,'average',names);

figure;
plot(tree);

%% Reroot tree w/ respect to 1991
selection = getbyname(tree,'99');

rootIndex = [];

for i = 1:length(selection)
    if selection(i) == 1
        rootIndex = i;
    end
end

rerootTree = reroot(tree,rootIndex);

figure;
plot(rerootTree)
 
%% Genetic distance vs. # of days
dist99 = 0;
dist00 = seqpdist([seq99(i99) seq00(i00)]);
dist01 = seqpdist([seq99(i99) seq01(i01)]);
dist03 = seqpdist([seq99(i99) seq03(i03)]);
dist06 = seqpdist([seq99(i99) seq06(i06)]);
dist07 = seqpdist([seq99(i99) seq07(i07)]);
dist08 = seqpdist([seq99(i99) seq08(i08)]);
dist10 = seqpdist([seq99(i99) seq10(i10)]);

xYears = [1999 2000 2001 2003 2006 2007 2008 2010];
yDist = [dist99 dist00 dist01 dist03 dist06 dist07 dist08 dist10];

figure;
scatter(xYears,yDist);
lsline;

figure;
scatter(xYears(2:end),yDist(2:end));
lsline;

%% Rate of evolution calculation
[eq erInf] = polyfit(xYears,yDist,1);
evRate = eq(1);
intercept = eq(2);

[eq2 erInf2] = polyfit(xYears(2:end),yDist(2:end),1);
evRate2 = eq2(1);
intercept2 = eq2(2);

%% Accuracy of evolution prediction
rand99 = randi(length(seq99)-1,1);
rand00 = randi(length(seq00)-1,1);
rand01 = randi(length(seq01)-1,1);
rand03 = randi(length(seq03)-1,1);
rand06 = randi(length(seq06)-1,1);
rand07 = randi(length(seq07)-1,1);
rand08 = randi(length(seq08)-1,1);
rand10 = randi(length(seq10)-1,1);

dt99 = seqpdist([seq99(i99) seq99(rand99)]);
dt00 = seqpdist([seq99(i99) seq00(rand00)]);
dt01 = seqpdist([seq99(i99) seq01(rand01)]);
dt03 = seqpdist([seq99(i99) seq03(rand03)]);
dt06 = seqpdist([seq99(i99) seq06(rand06)]);
dt07 = seqpdist([seq99(i99) seq07(rand07)]);
dt08 = seqpdist([seq99(i99) seq08(rand08)]);
dt10 = seqpdist([seq99(i99) seq10(rand10)]);

distSet = [dt99 dt00 dt01 dt03 dt06 dt07 dt08 dt10];

alpha = 0.05;

[distPred, error] = polyconf(eq2,xYears,erInf2,'alpha',alpha);

count = 0;

for i = 1:length(distSet)
    if distSet(i) <= distPred(i) + error(i) && distSet(i) >= distPred(i) - error(i)
        count = count + 1;
    end
end

accuracy = (count / length(distSet)) * 100;