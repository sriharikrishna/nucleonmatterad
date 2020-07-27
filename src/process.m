clear all

filename ='data_tap_all_dfo_snm';
n=2;
yscale='linear';

%filename ='data_tap_all_bfgs_snm';
%n=2;
%yscale='linear';

%filename ='data_tap_all_dfo_pnm';
%n=4;
%yscale='log';

%filename ='data_tap_all_bfgs_pnm';
%n=4;
%yscale='log';

dim=2+2*n;
fid = fopen([filename,'.txt']);

% Next section is a per-nucleus summary of the different Fayans computations
% that were run to collect all the required theoretical results.
j = 0;
while 1
    line = fgetl(fid);
    if line==-1
        break
    end
    b=strsplit(line);
    
    fprintf("length %d\n", length(b))
    if length(b)==2
        j=j+1;
        probstore(j,:)=cellfun(@str2num,b(1:end))';
        k=0;
    else
        k=k+1;
        datastore(k,:,j)=cellfun(@str2num,b(2:dim+1));
    end
end

for step=[0.1 0.2]
prob=0;
datastore(datastore==0)=NaN;
figure(1); clf; hold on;
for j=1:size(probstore,1)
    if probstore(j,2)==step
        prob=prob+1;
        figure(1)
        plot(datastore(:,2,j))
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
    end
end
figure(1);
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('f value')
title(sprintf('%s %f','DFO SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_1.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_1.png'])
end
#saveas(gcf, 'figure1.png')

figure(2)
plot(best)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('best f value')
title(sprintf('%s %f','DFO SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_2.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_2.png'])
end
#saveas(gcf, 'figure2.png')

figure(3)
plot(grad)
box on
axis tight
set(gca,'yscale','yscale')
xlabel('iterations')
ylabel('gradient norm')
title(sprintf('%s %f','DFO SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_3.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_3.png'])
end
#saveas(gcf, 'figure3.png')
end
clear all

%filename ='data_tap_all_dfo_snm';
%n=2;
%yscale='linear';

filename ='data_tap_all_bfgs_snm';
n=2;
yscale='linear';

%filename ='data_tap_all_dfo_pnm';
%n=4;
%yscale='log';

%filename ='data_tap_all_bfgs_pnm';
%n=4;
%yscale='log';

dim=2+2*n;
fid = fopen([filename,'.txt']);

% Next section is a per-nucleus summary of the different Fayans computations
% that were run to collect all the required theoretical results.
j = 0;
while 1
    line = fgetl(fid);
    if line==-1
        break
    end
    b=strsplit(line);
    
    fprintf("length %d\n", length(b))
    if length(b)==2
        j=j+1;
        probstore(j,:)=cellfun(@str2num,b(1:end))';
        k=0;
    else
        k=k+1;
        datastore(k,:,j)=cellfun(@str2num,b(2:dim+1));
    end
end

for step=[0.1 0.2]
prob=0;
datastore(datastore==0)=NaN;
figure(1); clf; hold on;


for j=1:size(probstore,1)
    if probstore(j,2)==step
        prob=prob+1;
        figure(1)
#plot(3:20, datastore(3:end,2,j))
        plot(datastore(:,2,j))
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
    end
end
figure(1);
box on
axis tight
#set(gca,'yscale','log')
#set(gca, 'YTick', [10.^-2 10.^-1 10.^0 10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8])
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('f value')
#title(sprintf('%s %f','BFGS SNM (ITERATION 3 ONWARDS)', step))
title(sprintf('%s %f','BFGS SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_1.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_1.png'])
end
#saveas(gcf, 'figure1.png')

figure(2)
plot(best)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('best f value')
title(sprintf('%s %f','BFGS SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_2.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_2.png'])
end
#saveas(gcf, 'figure2.png')

figure(3)
plot(grad)
box on
axis tight
set(gca,'yscale','log')
set(gca, 'YTick', [10.^-2 10.^-1 10.^0 10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8])
xlabel('iterations')
ylabel('gradient norm')
title(sprintf('%s %f','BFGS SNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_3.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_3.png'])
end
#saveas(gcf, 'figure3.png')
end
clear all

%filename ='data_tap_all_dfo_snm';
%n=2;
%yscale='linear';

%filename ='data_tap_all_bfgs_snm';
%n=2;
%yscale='linear';

filename ='data_tap_all_dfo_pnm';
n=4;
yscale='log';

%filename ='data_tap_all_bfgs_pnm';
%n=4;
%yscale='log';

dim=2+2*n;
fid = fopen([filename,'.txt']);

% Next section is a per-nucleus summary of the different Fayans computations
% that were run to collect all the required theoretical results.
j = 0;
while 1
    line = fgetl(fid);
    if line==-1
        break
    end
    b=strsplit(line);
    
    fprintf("length %d\n", length(b))
    if length(b)==2
        j=j+1;
        probstore(j,:)=cellfun(@str2num,b(1:end))';
        k=0;
    else
        k=k+1;
        datastore(k,:,j)=cellfun(@str2num,b(2:dim+1));
    end
end

for step=[0.1 0.2]
prob=0;
datastore(datastore==0)=NaN;
figure(1); clf; hold on;
for j=1:size(probstore,1)
    if probstore(j,2)==step
        prob=prob+1;
        figure(1)
        plot(datastore(:,2,j))
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
    end
end
figure(1);
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('f value')
title(sprintf('%s %f','DFO PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_1.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_1.png'])
end
#saveas(gcf, 'figure1.png')

figure(2)
plot(best)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('best f value')
title(sprintf('%s %f','DFO PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_2.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_2.png'])
end
#saveas(gcf, 'figure2.png')

figure(3)
plot(grad)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('gradient norm')
title(sprintf('%s %f','DFO PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_3.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_3.png'])
end
#saveas(gcf, 'figure3.png')
end
clear all

%filename ='data_tap_all_dfo_snm';
%n=2;
%yscale='linear';

%filename ='data_tap_all_bfgs_snm';
%n=2;
%yscale='linear';

%filename ='data_tap_all_dfo_pnm';
%n=4;
%yscale='log';

filename ='data_tap_all_bfgs_pnm';
n=4;
yscale='log';

dim=2+2*n;
fid = fopen([filename,'.txt']);

% Next section is a per-nucleus summary of the different Fayans computations
% that were run to collect all the required theoretical results.
j = 0;
while 1
    line = fgetl(fid);
    if line==-1
        break
    end
    b=strsplit(line);
    
    fprintf("length %d\n", length(b))
    if length(b)==2
        j=j+1;
        probstore(j,:)=cellfun(@str2num,b(1:end))';
        k=0;
    else
        k=k+1;
        datastore(k,:,j)=cellfun(@str2num,b(2:dim+1));
    end
end

for step=[0.1 0.2]
prob=0;
datastore(datastore==0)=NaN;
figure(1); clf; hold on;
for j=1:size(probstore,1)
    if probstore(j,2)==step
        prob=prob+1;
        figure(1)
        plot(datastore(:,2,j))
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
    end
end
figure(1);
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('f value')
title(sprintf('%s %f','BFGS PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_1.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_1.png'])
end
#saveas(gcf, 'figure1.png')

figure(2)
plot(best)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('best f value')
title(sprintf('%s %f','BFGS PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_2.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_2.png'])
end
#saveas(gcf, 'figure2.png')

figure(3)
plot(grad)
box on
axis tight
set(gca,'yscale',yscale)
xlabel('iterations')
ylabel('gradient norm')
title(sprintf('%s %f','BFGS PNM ', step))
if step==0.1
print('-dpng',['fig_',filename,'_','0_1','_3.png'])
end
if step==0.2
print('-dpng',['fig_',filename,'_','0_2','_3.png'])
end
#saveas(gcf, 'figure3.png')
end

