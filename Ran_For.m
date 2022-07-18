% clc
% clear all
function predict_label = Ran_For(FV,D)

if D == 10
     load('train_data1_10D.mat');
     train_data = train_data1_10D;
elseif D == 30
     load('train_data1_30D.mat');
     train_data = train_data1_30D;
elseif D == 50
    load('train_data1_50D.mat');
    train_data = train_data1_50D;
elseif D == 100
    load('train_data1_100D.mat');
    train_data = train_data1_100D;
end

nTree = 10;
% Label = randi(3,1,30);
% Label = [1 1 1 3 3 1 1 2 1 3 2 1 3 1 1 3 2 2 2 3 2 1 1 1 1 1 1 1 2 3];
% Feature_Label = HE;
% Feature_Label = [HE FDC Label'];
train_label = train_data(:,10);
train_data = train_data(:,1:9);

test_data = FV;

% train_data = Feature_Label(1:15,1:10);
% test_data = Feature_Label(16:end,1:10);
% test_data = train_data(randperm(30,10),:);
% test_data(randperm(100,10)) = unidrnd(3);
% train_label = Feature_Label(1:15,11);
B = TreeBagger(nTree,train_data,train_label,'method','classification');
predict_label = predict(B,test_data);
