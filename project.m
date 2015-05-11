%% Load in sequences data file
filename = fopen('/Users/atsiros/seqData.csv');
cols = textscan(filename,'%f %s %s %s %s','delimiter',',');
fclose(filename);

%% Define columns
sequences = cols{5};
sequences(456) = [];
type = cols{4};
type(456) = [];
date = cols{2};
time = cols{1};

%% Set up training data
seqTrain = sequences(2:410);
seqTrain(370) = [];
seqTrain(368) = [];
seqTrain(362) = [];
seqTrain(357) = [];
seqTrain(225) = [];
seqTrain(177) = [];
seqTrain(176) = [];
seqTrain(168) = [];
seqTrain(167) = [];
seqTrain(166) = [];
seqTrain(163) = [];
seqTrain(41) = [];
seqTrain(8) = [];

typeTrain = type(2:410);
typeTrain(370) = [];
typeTrain(368) = [];
typeTrain(362) = [];
typeTrain(357) = [];
typeTrain(225) = [];
typeTrain(177) = [];
typeTrain(176) = [];
typeTrain(168) = [];
typeTrain(167) = [];
typeTrain(166) = [];
typeTrain(163) = [];
typeTrain(41) = [];
typeTrain(8) = [];

%% Set up testing data
seqTest = sequences(397:end);
typeTest = type(397:end);

%% Create multiple sequence alignment and define length
MSA = multialign(seqTrain);

l = length(MSA);

%% Build PHMM model
model = hmmprofstruct(l);
model2 = hmmprofestimate(model,MSA);

%% Training
A = [];
B = [];
C = [];

for x = 1:length(typeTrain)
    if typeTrain{x} == 'a'
        score = hmmprofalign(model2,seqTrain{x});
        A = [A; score];
    elseif typeTrain{x} == 'b'
        score = hmmprofalign(model2,seqTrain{x});
        B = [B; score];
    elseif typeTrain{x} == 'c'
        score = hmmprofalign(model2,seqTrain{x});
        C = [C; score];
    else
        score = hmmprofalign(model2,seqTrain{x});
        G = [G; score];
    end
end

%% 95% Confidence intervals for a & b subtypes
ciB = [];
upperB = mean(B) + std(B)*1.96;
lowerB = mean(B) - std(B)*1.96;
ciB = [lowerB upperB];

ciA = [];
upperA = mean(A) + std(A)*1.96;
lowerA = mean(A) - std(A)*1.96;
ciA = [lowerA upperA];

%% Testing
testResultsA = [];
testResultsB = [];
testAssignments = [];

for x = 1:length(seqTest)
    score = hmmprofalign(model2,seqTest{x});
    if score >= ciA(1) && score <= ciA(2)
        testAssignments = [testAssignments; 'a'];
        testResultsA = [testResultsA; score];
    elseif score >= ciB(1) && score <= ciB(2)
        testAssignments = [testAssignments; 'b'];
        testResultsB = [testResultsB; score];
    else
        testAssignments = [testAssignments; 'x'];
    end
end

%% Results
count = 0;

for x = 1:length(typeTest)
    if typeTest{x} == testAssignments(x)
        count = count + 1;
    end
end

accuracy = count / length(typeTest);

tp = 0;
fp = 0;
fn = 0;

for x = 1:length(typeTest)
    if typeTest{x} == testAssignments(x) && typeTest{x} == 'b'
        tp = tp + 1;
    elseif typeTest{x} ~= testAssignments(x) && typeTest{x} == 'b'
        fp = fp + 1;
    elseif typeTest{x} ~= testAssignments(x) && typeTest{x} == 'a'
        fn = fn + 1;
    end
end

precision = tp / (tp + fp);
recall = tp / (tp + fn);