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

% figure; scatter(MidPoint_Date_UTC,taxonabundance)
%plots sample abundances against time for visual check

firstyear=1990;
%thru
lastyear=2020;
%example period

southlim=35;
northlim=65;
westlim=75; %in degrees west of greenwich
eastlim=-23; %in degrees east of greenwich
%total extent of CPR BEAMS data

useindex=find((lat>=southlim).*(lat<northlim).*(long>=-westlim).*(long<eastlim).*(year>=firstyear).*(year<=lastyear));
%first restricts the data to a box and a time period

polylong=[-71.4 -70.4 -70.4 -71.4];
polylat=[39 39.5 41.5 41.5];
%NES LTER region

in = inpolygon(long(useindex),lat(useindex),polylong,polylat);
% in = inpolygon(xq,yq,xv,yv) returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
% note that the polygon is drawn on a flat plate projection, not on the
% globe, and the edges do not correspond to the shortest distance between
% vertices on the globe.
% [in,on] = inpolygon(xq,yq,xv,yv) also returns on indicating if the query points are on the edge of the polygon area.

% figure; line(polylong,polylat)
% hold on; scatter(long(useindex(in)),lat(useindex(in)))
% 
% figure; scatter(MidPoint_Date_UTC(useindex(in)),taxonabundance(useindex(in)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sets up data samples in a compact cell array

yearindex=year(useindex(in))+1-firstyear;
monthindex=month(useindex(in));
%transforms the sample time and location data into reference indices

abund=taxonabundance(useindex(in));

for y=1:lastyear+1-firstyear
    for m=1:12
                dataset{y,m}=[];
                datasetinfo{y,m}=[firstyear+y-1 m];
    end
end

for n=1:length(yearindex)
    dataset{yearindex(n),monthindex(n)}=[dataset{yearindex(n),monthindex(n)} abund(n)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%monthly timeseries with NaNs

for y=1:lastyear+1-firstyear
    for m=1:12
                abundseries((y-1)*12+m)=nanmean(dataset{y,m}); %monthly timeseries with NaNs
                seriesyear((y-1)*12+m)=datasetinfo{y,m}(1);
                seriesmonth((y-1)*12+m)=datasetinfo{y,m}(2);
    end
end

%monthly medians
for m=1:12
            mabundseries(m)=nanmedian(abundseries(([1:lastyear+1-firstyear]-1)*12+m)); %monthly medians
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fill timeseries

for y=1:lastyear+1-firstyear
    for m=1:12
                if isnan(abundseries((y-1)*12+m))
                    nnabundseries((y-1)*12+m)=mabundseries(m);     %fill timeseries with median for that location and month of year

                else
                    nnabundseries((y-1)*12+m)=abundseries((y-1)*12+m); %use actual monthly data where available
                end
            
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
subplot(2,2,1); %a map of the defined area and samples
hold on; 

line([polylong, polylong(1)],[polylat, polylat(1)],'Color','m')

scatter(long(useindex(in)),lat(useindex(in)))

load coastlines
% axes('position',[0 0 1 1]); hold on
[latcells, loncells] = polysplit(coastlat, coastlon);
for n=1:length(latcells)
    fill(loncells{n},latcells{n},'g');
end
for n=[19 20 21 22 23 95 100 101 106 107]
    fill(loncells{n},latcells{n},[1 1 1]);
end
%add a coastline for reference

xlabel('longitude')
ylabel('latitiude')
xlim([min(polylong)-1 max(polylong)+1])
ylim([min(polylat)-1 max(polylat)+1])
title('samples')

subplot(2,2,2) %monthly abundances constructed from samples in the area
            plot(firstyear+1/12:1/12:lastyear+1,abundseries); 
xlabel('year')
ylabel('abundance')
title(taxon_selection)

subplot(2,2,3) %median monthly values
            plot(mabundseries,'r');
xlabel('month')
ylabel('median abundance')
title('seasonal profile')         
            
subplot(2,2,4) %filled time-series
hold on;
            plot(firstyear+1/12:1/12:lastyear+1,nnabundseries,'r');
            plot(firstyear+1/12:1/12:lastyear+1,abundseries,'b.-'); 
            xlim([firstyear-1 lastyear+1])
 xlabel('year')
ylabel('abundance (missing values filled)')
title('time-series')

set(gcf, 'paperpositionmode','manual','paperunits','inches','paperposition',[0 0 10 10],'papersize',[10 10])
print(gcf,'-djpeg', '-r300', 'CPRBEAMStimeseriesexample_arearestriction_2022.jpg')
% close all

% %these variables can be used for further analysis

% dataset
% %lists the abundances found in each month examined

% datasetinfo
% %records the year and month associated with each of the abundances in
% %dataset

% abundseries
% %mean of the abundances found in each month examined, NaN where
% %not available

% mabundseries
% %median of the mean monthly abundances found (a seasonal profile), NaN where
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
