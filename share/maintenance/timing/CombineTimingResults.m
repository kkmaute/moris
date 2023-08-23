function main

clear all
close all

titan              =''; %'/home/maute/codes/moris_titan/build_opt';
titan_xtk_refactor =''; %'/home/maute/codes/moris_titan_xtk_refactor/build_opt';
github             = [getenv('MORISROOT'), '/build_opt/TimingResults'];

cd(github)
if not(isfolder('plots'))
    mkdir('plots')
end

testcases=ls('*.timing');
testcases=split(testcases);

for i=1:length(testcases)
    
    if size(testcases{i})==0; continue;end
    
    num_titan=0;
    mat_titan=[];
    if length(titan)>0
        fname=sprintf('%s/%s',titan,testcases{i});
        mat_titan=textread(fname);
        num_titan=size(mat_titan,1);
    end
    
    num_titan_xtk_refactor=0;
    mat_titan_xtk_refactor=[];
    if length(titan_xtk_refactor) > 0
        fname=sprintf('%s/%s',titan_xtk_refactor,testcases{i});
        mat_titan_xtk_refactor=textread(fname);
        num_titan_xtk_refactor=size(mat_titan_xtk_refactor,1);
    end
    
    fname=sprintf('%s/%s',github,testcases{i});
    mat_github=textread(fname);
    num_github=size(mat_github,1);
    
    timings=[mat_titan;mat_titan_xtk_refactor;mat_github];
    
    figure(1)
    clf
    plot(timings(:,2)); hold on
    plot( [num_titan num_titan],[0 max(timings(:,2))],'k-'); hold on
    plot( [num_titan+num_titan_xtk_refactor num_titan+num_titan_xtk_refactor],[0 max(timings(:,2))],'k-');
    title(replace(testcases{i},'_','-'));
    saveas(gcf, ['plots/', testcases{i}, '_F1.png']);
    
    figure(2)
    clf
    if length(titan)>0
        plotdays(mat_titan,'ks'); hold on
    end
    if length(titan_xtk_refactor) >0
        plotdays(mat_titan_xtk_refactor,'b*'); hold on
    end
    plotdays(mat_github,'rs'); hold on
    title(replace(testcases{i},'_','-'));
    saveas(gcf, ['plots/', testcases{i}, '_F2.png']);
end
end

function plotdays(mat,str)

days=convdate(mat(:,1));
[days,itx]=sort(days);
tims=mat(itx,2);
plot(days',tims,str); hold on

end

function [outvec]=convdate(datevec)

for i=1:length(datevec)
    
    digits=dec2base(datevec(i),10) - '0';
    
    outvec(i) = datenum(datetime(2000+10*digits(1)+digits(2),10*digits(3)+digits(4),10*digits(5)+digits(6)))-datenum(datetime(2020,1,1));
    
end


end