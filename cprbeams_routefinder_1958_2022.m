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

% routechoice='Z-'
% routechoice='ZB'
% routechoice='ZC'
% routechoice='MB'
routechoice='MC'
% routechoice='MD'
% routechoice='EB'
% routechoice='NWP'
% routechoice='MC'
% routechoice='BC'
% routechoice='BD'

routecheck=strfind(SampleId,routechoice);
routeidlist=[];
for n=1:length(routecheck)
    if length(routecheck{n})>0
        routeidlist=[routeidlist n];
    end
end
%routeidlist provides a list of indices of all samples in data that belong
%to the route choice
routesamples=SampleId(routeidlist);
routelats=lat(routeidlist);
routelongs=long(routeidlist);
routedates=MidPoint_Date_UTC(routeidlist);

figure; 
subplot(1,2,1); scatter(routelongs,routelats); title(routechoice); xlabel('Longitude'); ylabel('Latitude')
subplot(1,2,2); plot(sort(routedates),1:length(routedates)); xlabel('Date'); ylabel('Accumulated samples')