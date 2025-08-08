close all
clear all

%data released in 2025 is current on BCO-DMO: this is all data processed
%by end 2024, covering 1958-2022 inclusive
metadata=readtable('https://datadocs.bco-dmo.org/file/DwDK8ZqH3LYvDy/765141_v6_cpr-list-taxa.csv');
data=readtable('https://datadocs.bco-dmo.org/file/WWrAqKPH6qLEvw/765141_v6_cpr-plankton-abundance.csv');

accepted_id=metadata.Accepted_ID;
%contains a list of all taxa identities according to CPR numbering
%convention
aphia_id=metadata.Aphia_ID;
%contains a list of all taxa identities according to aphia numbering
%convention
name_cpr=metadata.Taxon_Name;
%contains a list of all taxa identities according to cpr naming convention
name_worms=metadata.WoRMS_Name;
%contains a list of all taxa identities according to worms naming convention
DRI=metadata.DRI;
%contains a list of first inclusion dates (taxon would not have been looked
%for before this date)
counting_method=metadata.counting_method;
%1 for large zooplankton (eye count), 2 for small zoomplankton (semi-quantitative traverse count), 3 for phytoplankton (semi-quantitative traverse count)

SampleId=data.SampleId;
lat=data.Latitude;
long=data.Longitude;
MidPoint_Date_UTC=datetime(data.MidPoint_Date_UTC,'format','yyyy-MM-dd''T''HH:mm''Z');
%sample id, postion and time (GMT)

year=data.Year;
month=data.Month;
day=data.Day;
hour=data.Hour;
%redundant time information

%EXAMPLE
worms_taxon_selection='Oithona';
name_cpr(find(strcmp(name_worms,worms_taxon_selection)))
%identifies CPR taxon names associated with this worms taxon name (there may be
%more than one)

taxon_selection='Oithona spp.';
%selected CPR taxon name
selectedid=['id_', num2str(metadata.Accepted_ID(find(strcmp(name_cpr,taxon_selection))))];
%finds the CPR id associated with this CPR taxon name

eval(['taxonabundance=data.' selectedid ';']);
%extracts abundance

figure; scatter(MidPoint_Date_UTC,taxonabundance)
%plots sample abundances against time for visual check

firstyear=1990;
lastyear=2022;
southlim=35;
northlim=65;
westlim=75; %in degrees west
eastlim=-23; %in degrees east

useindex=find((lat>=southlim).*(lat<northlim).*(long>=-westlim).*(long<eastlim).*(year>=firstyear).*(year<=lastyear));
%restricts the data to a box


squaredeg=2;
nanthreshold=0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sets up data samples in a compact cell array

yearindex=year(useindex)+1-firstyear;
monthindex=month(useindex);
latindex=floor((lat(useindex)-southlim)/squaredeg)+1;
longindex=floor((long(useindex)+westlim)/squaredeg)+1;
%transforms the sample time and location data into reference indices

abund=taxonabundance(useindex);

lat1max=max(latindex);
long1max=max(longindex);
for y=1:lastyear+1-firstyear
    for m=1:12
        for lat1=1:lat1max
            for long1=1:long1max
                dataset{y,m,lat1,long1}=[];
                datasetinfo{y,m,lat1,long1}=[firstyear+y-1 m southlim+(lat1-1)*squaredeg -westlim+(long1-1)*squaredeg];
            end
        end
    end
end

for n=1:length(yearindex)
    dataset{yearindex(n),monthindex(n),latindex(n),longindex(n)}=[dataset{yearindex(n),monthindex(n),latindex(n),longindex(n)} abund(n)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%monthly timeseries with NaNs
for y=1:lastyear+1-firstyear
    for m=1:12
        for lat1=1:lat1max
            for long1=1:long1max
                abundseries{lat1,long1}((y-1)*12+m)=nanmean(dataset{y,m,lat1,long1}); %monthly timeseries with NaNs
                seriesyear{lat1,long1}((y-1)*12+m)=datasetinfo{y,m,lat1,long1}(1);
                seriesmonth{lat1,long1}((y-1)*12+m)=datasetinfo{y,m,lat1,long1}(2);
                serieslat{lat1,long1}((y-1)*12+m)=datasetinfo{y,m,lat1,long1}(3);
                serieslong{lat1,long1}((y-1)*12+m)=datasetinfo{y,m,lat1,long1}(4);
            end
        end
    end
end
figure; hold on;
load coastlines
axes('position',[0 0 1 1]); hold on
[latcells, loncells] = polysplit(coastlat, coastlon);
for n=1:length(latcells)
    fill(loncells{n},latcells{n},'g');
end
for n=[19 20 21 22 23 95 100 101 106 107]
    fill(loncells{n},latcells{n},[1 1 1]);
end
xlabel('longitude')
ylabel('latitiude')
xlim([-westlim eastlim])
ylim([southlim northlim])
scatter(long(useindex),lat(useindex),'r.')
for lat1=1:lat1max
    for long1=1:long1max
        if sum(isnan(abundseries{lat1,long1}))<nanthreshold*length(abundseries{lat1,long1});
            axes('position',[(((long1-1)*squaredeg-westlim)+westlim)/(eastlim+westlim) (((lat1-1)*squaredeg+southlim)-southlim)/(northlim-southlim) (squaredeg-0.5)/(eastlim+westlim) (squaredeg-0.5)/(northlim-southlim)])
            hold on;
            plot(firstyear+1/12:1/12:lastyear+1,abundseries{lat1,long1}); title([num2str((lat1-1)*squaredeg+southlim) ' ' num2str((long1-1)*squaredeg-westlim)])
        end
    end
end
set(gcf, 'paperpositionmode','manual','paperunits','inches','paperposition',[0 0 eastlim+westlim northlim-southlim],'papersize',[eastlim+westlim northlim-southlim])
print(gcf,'-djpeg', '-r300', 'CPRBEAMSdataexample_2022.jpg')

%monthly medians
for m=1:12
    for lat1=1:lat1max
        for long1=1:long1max
            mabundseries{lat1,long1}(m)=nanmedian(abundseries{lat1,long1}(([1:lastyear+1-firstyear]-1)*12+m)); %monthly medians

        end
    end
end
figure; hold on;
load coastlines
axes('position',[0 0 1 1]); hold on
[latcells, loncells] = polysplit(coastlat, coastlon);
for n=1:length(latcells)
    fill(loncells{n},latcells{n},'g');
end
for n=[19 20 21 22 23 95 100 101 106 107]
    fill(loncells{n},latcells{n},[1 1 1]);
end
xlabel('longitude')
ylabel('latitiude')
xlim([-westlim eastlim])
ylim([southlim northlim])
for lat1=1:lat1max
    for long1=1:long1max
        if sum(isnan(mabundseries{lat1,long1}))<nanthreshold*length(mabundseries{lat1,long1});
            hold on;
            axes('position',[(((long1-1)*squaredeg-westlim)+westlim)/(eastlim+westlim) (((lat1-1)*squaredeg+southlim)-southlim)/(northlim-southlim) (squaredeg-0.5)/(eastlim+westlim) (squaredeg-0.5)/(northlim-southlim)])
            hold on;
            plot(mabundseries{lat1,long1},'r');
            
        end
    end
end
set(gcf, 'paperpositionmode','manual','paperunits','inches','paperposition',[0 0 eastlim+westlim northlim-southlim],'papersize',[eastlim+westlim northlim-southlim])
print(gcf,'-djpeg', '-r300', 'CPRBEAMSmonthlymediansexample_2022.jpg')

%     fill timeseries
for y=1:lastyear+1-firstyear
    for m=1:12
        for lat1=1:lat1max
            for long1=1:long1max
                if isnan(abundseries{lat1,long1}((y-1)*12+m))
                    nnabundseries{lat1,long1}((y-1)*12+m)=mabundseries{lat1,long1}(m);     %fill timeseries with median for that location and month of year

                else
                    nnabundseries{lat1,long1}((y-1)*12+m)=abundseries{lat1,long1}((y-1)*12+m); %use actual monthly data where available
                end
            end
        end
    end
end

figure; hold on;
load coastlines
axes('position',[0 0 1 1]); hold on
[latcells, loncells] = polysplit(coastlat, coastlon);
for n=1:length(latcells)
    fill(loncells{n},latcells{n},'g');
end
for n=[19 20 21 22 23 95 100 101 106 107]
    fill(loncells{n},latcells{n},[1 1 1]);
end
xlabel('longitude')
ylabel('latitiude')
xlim([-westlim eastlim])
ylim([southlim northlim])
for lat1=1:lat1max
    for long1=1:long1max
        if sum(isnan(abundseries{lat1,long1}))<nanthreshold*length(abundseries{lat1,long1});
            axes('position',[(((long1-1)*squaredeg-westlim)+westlim)/(eastlim+westlim) (((lat1-1)*squaredeg+southlim)-southlim)/(northlim-southlim) (squaredeg-0.5)/(eastlim+westlim) (squaredeg-0.5)/(northlim-southlim)])
            hold on;
            plot(firstyear+1/12:1/12:lastyear+1,nnabundseries{lat1,long1},'r');
            plot(firstyear+1/12:1/12:lastyear+1,abundseries{lat1,long1},'b.-'); title([num2str((lat1-1)*squaredeg+southlim) ' ' num2str((long1-1)*squaredeg-westlim)])
            xlim([firstyear-1 lastyear+1])
        end
    end
end
set(gcf, 'paperpositionmode','manual','paperunits','inches','paperposition',[0 0 eastlim+westlim northlim-southlim],'papersize',[eastlim+westlim northlim-southlim])
print(gcf,'-djpeg', '-r300', 'CPRBEAMStimeseriesexample_2022.jpg')
% close all


% %these variables can be used for further analysis

% dataset
% %lists the abundances found in each year/month/lat/long examined

% datasetinfo
% %records the year, month, lat and long (bottom left corner of square from which the data is taken) associated with each of the abundances in
% %dataset

% abundseries
% %mean of the abundances found in each square in each month examined, NaN where
% %not available

% mabundseries
% %median of the mean monthly abundances found (a seasonal profile) in each square, NaN where
% %not available

% nnabundseries
% %time-series made by filling missing values in abundseries with
% %mabundseries values where available

% seriesyear
% %records the year associated with each of the abundances in
% % abundseries and nnabundseries

% seriesmonth
% %records the month associated with each of the abundances in
% % abundseries and nnabundseries

% serieslat
% %records the lat of the bottom left corner of the square associated with each of the abundances in
% % abundseries and nnabundseries

% serieslong
% %records the long of the bottom left corner of the square associated with each of the abundances in
% % abundseries and nnabundseries

