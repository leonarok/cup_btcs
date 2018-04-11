clear all
close all
clc

npi=52;
npj=102;

tend=3600;
dt=0.15;
printTimes=100;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);

fileloc = 'output/temp/temp_    .dat';
count = '    ';

tempMid=zeros(1,printSteps);tempTopRight=tempMid;
tempMean=tempMid;tempMidTop=tempMid;


time=print_dt:print_dt:printSteps*print_dt;

for n=1:printSteps
    
    fileTime=num2str(n*print_dt);
    
    if length(fileTime)==2
        count(3:4)=fileTime;
    elseif length(fileTime)==3
        count(2:4)=fileTime;
    elseif length(fileTime)==4
        count(1:4)=fileTime;
    end
    fileloc(18:21)=count;
    a=dlmread(fileloc);
    
    tempMid(n)=a(npi/2,npj/2);
    tempMidTop(n)=a(npi/2,ceil(npj*3/4));
    tempMean(n)=mean(mean(a(2:npi-1,2:npj-1)));
    tempTopRight(n)=a(ceil(npi*3/4),ceil(npj*3/4));
    
end

hold on
plot(time,tempMean)
plot(time,tempMid)
plot(time,tempMidTop)
plot(time,tempTopRight)
legend('mean temp','temp at midnode','mid top temp','top right temp')
