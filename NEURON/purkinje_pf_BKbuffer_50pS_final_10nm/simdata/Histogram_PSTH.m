%% 11 September 2024 %% Gabriela Cirtala
% This code will load the data resulting from Yunliang's paper
% from 500 simulations (Vsoma) and will plot the PSTH

clear all
close all

%% Define branch and corresponding number of PF
br = 12;
nbPF = 0; %fake number of PF in Yunliang's code
tnbPF = nbPF*5;%actual number of PF in Yunliang's code is multiplied by 5
%% Define time in ms:
t=0:0.02:700;

pathnumb = ['fig5_50pS_final_10nm_-0.1316inj']
cd (pathnumb)


%% define bin width for the histogram: 
binsize = 2;

%% Define number of trials:
iter_list =1:500;
v1vec = []; spikelist = [];

%% Load trials to obtain the voltage trails (at soma)
 sn=['003_' num2str(br,'%03.f') '_' num2str(br,'%03.f') '_' num2str(nbPF,'%03.f') '_'];
 for iter = iter_list
          v1=load([sn num2str(iter,'%03.f') '_vsoma.dat']);
          v1vec=[v1vec v1];
 end

 %% Finds spikes and makes sure that they are not multiply defined
for j=1:500
x=find(v1vec(:,j)>-10);
for i=1:(length(x)-1)
if x(i+1)-x(i)<2
    x(i+1)=0;
end
end
x(x==0)=[];
for i=1:(length(x)-1)
if x(i+1)-x(i)<10
   x(i+1)=0;
end
end
x(x==0)=[];

for i=1:(length(x)-1)
if x(i+1)-x(i)<10
   x(i+1)=0;
end
end
x(x==0)=[];
st = x;
spikelist = [spikelist;x];
end

tspikelist =t(spikelist);

 %% draw histogram with binwidth specified above
 figure(1)
 h = histogram(tspikelist,'BinWidth',binsize,'FaceColor','k');
 x_bin = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
 values_bin = h.Values; %% contains the values corresponding to the bin
 set(gca,'FontSize',25)
 xlabel('time (ms)')
 ylabel('Spiking rate')
 xlim([50,700])
 ylim([0,300])
 %%hgexport(gcf, 'PSTH_50ps', hgexport('factorystyle'), 'Format', 'png','Resolution','500')

%% save results
save spikes.mat tspikelist
save values_histogram.mat values_bin

disp(mean(tspikelist))

cd ('..')

