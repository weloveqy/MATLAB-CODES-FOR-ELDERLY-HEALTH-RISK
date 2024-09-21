%% codes for article entitled 'Intensifying Health Risk Inequality for the Elderly'
%% author: Yuan Qi  ; Email : qqyy@hku.hk
%% Step 1:calculate the ATEI; Step 2: calculate the CHRI ；  Step 3: tests


%% Step 1:calculate the ATEI ....

%% heat waves and cold spells in the past 31 years and future 30 years
% read and input
clear

listing=dir(' \ERA5 2m air temperature\*.mat'); list_size=size(listing,1);
sort_nat_name=sort_nat({listing.name}); % sort_nat: Natural Order Sort
cycle_count=0;
for l=1:list_size
    cycle_count=cycle_count+1; cycle_count
    name=sort_nat_name(l); name=char(name); filepath=[' \ERA5 2m air temperature\' name];
    Data=load(filepath); Tmax=Data.Tmax_matrix; 
    Tmax_Everyday(:,:,l)=single(Tmax); 
end

 data_01=load(' \Auxiliary data\Land_Water_01_Global_Grid_0.25deg.mat'); land_water=data_01.data; 
 land_water=flipud(land_water); land_water=rot90(land_water,3);
 winter_threshold_1=load(' \Winter_Threshold_32.mat'); winter_threshold=winter_threshold_1.winter_threshold;
 summer_threshold_1=load(' \Summer_Threshold_33.mat'); summer_threshold=summer_threshold_1.summer_threshold;
 winter_threshold=single(winter_threshold); summer_threshold=single(summer_threshold);

% extreme temperature events definition, intensity and space constrains of T
% start_day=date1; duration=152; Year_num=33;   % month 5 6 7 8 9; 11 12 1 2 3; 1990-2022  33/32
 start_day=date2; duration=151; Year_num=32;    % select one of the start_day & duration
 HW_central=single(zeros(1439, 720,duration, Year_num)); HW_neighbor=single(zeros(1439, 720,duration, Year_num)); % heat waves
 CW_central=single(zeros(1439, 720,duration, Year_num)); CW_neighbor=single(zeros(1439, 720,duration, Year_num)); % cold spells
 
 cycle_count=0; 
for i=2:1439
     for j=2:720  % north hemisphere hot/cold waves 2-453
        day_count=0; cycle_count=cycle_count+1;
        cycle_count  
        if land_water(i,j-1)==1   
            for d=start_day:start_day+duration-1  %calendar day May-August 121-243; November-February
                cout_num=0; day_count=day_count+1; Sliding_cell_total=[];
        
                for m=1:Year_num   % total 33 years             
                  for n=d-3+365*(m-1):d+3+365*(m-1)  % total 7 days  sliding window     
                     cout_num=cout_num+1;
                     Neighbor_grid_cells=Tmax_Everyday(i-1:i+1,j-1:j+1,n); 
                     Sliding_cell_total(:,cout_num)=Neighbor_grid_cells(:);  % 33*7*9=2079
                  end        
                end    
            
                Threshold_h_90th=prctile(Sliding_cell_total(:),90); % Threshold_h_95th=prctile(Sliding_cell_total(:),95); % select one of the Threshold groups
                Threshold_c_10th=prctile(Sliding_cell_total(:),10); % Threshold_c_05th=prctile(Sliding_cell_total(:),05);
        
                for mm=1:Year_num              
                    Central_grid=Tmax_Everyday(i,j,d+365*(mm-1));                
                    Neighbor_grid_cells=[];  Neighbor_grid_cells=Tmax_Everyday(i-1:i+1,j-1:j+1,d+365*(mm-1));
                    threshold_sum=summer_threshold(i,j,mm);  threshold_win=winter_threshold(i,j,mm);
                    if Central_grid > threshold_sum                    
                        Intensity=abs(Central_grid/Threshold_h_90th); HWs_count=length(Neighbor_grid_cells>Threshold_h_90th); %  heat waves
                        HW_central(i, j, day_count, mm)=Intensity; HW_neighbor(i, j, day_count, mm)=HWs_count;   
                    end
                    if Central_grid < threshold_win     
                        Intensity=abs(Threshold_c_10th/Central_grid); CWs_count=length(Neighbor_grid_cells<Threshold_c_10th); %  cold spells
                        CW_central(i, j, day_count, mm)=Intensity; CW_neighbor(i, j, day_count, mm)=CWs_count;
                    end
                end   
         
            end
        else
            continue % loop to next grid
        end
    end
end

%% if 3 concective days, time constrains of T, get acumulative intensity of every event
Intensity_grid=zeros(1439, 453,Year_num); cycle_count=0; Intensity_counts=zeros(1439, 453,Year_num);
for i=2:1439
    for j=2:720
        cycle_count=cycle_count+1; cycle_count
        for k=1:Year_num
           %  Single_year_series= HW_central(i,j,:,k);   Single_year_neibour= HW_neighbor(i,j,:,k);
            Single_year_series= CW_central(i,j,:,k);   Single_year_neibour= CW_neighbor(i,j,:,k);
            Single_year_series=Single_year_series(:);  Single_year_neibour= Single_year_neibour(:);
            Cum_Int=0; Cum_Int_2=0; count_l=0; Cum_counts=0; Cum_counts_2=0;
            for l=1:length(Single_year_series)
                if Single_year_series(l)>1 && Single_year_neibour(l)>1  %  hws cws 
                    count_l=count_l+1; 
                    Cum_Int=Cum_Int+1;  Cum_counts=Cum_counts+Single_year_series(l); 
                    if count_l>2        
                       Cum_Int_2=Cum_Int_2+Cum_Int; Cum_Int=0;
                       Cum_counts_2=Cum_counts_2+Cum_counts; Cum_counts=0;
                    end
                else
                    count_l=0; Cum_Int=0;  Cum_counts=0;
                end     
            end
            Intensity_grid(i,j,k)=Cum_Int_2; Intensity_counts(i,j,k)=Cum_counts_2;
        end
    end
end

% calculate total intensity of 31yrs
Intensity_value=zeros(1439,720); Counts_value=zeros(1439,720);
for i=2:1439
    for j=2:720
        Grid_value=Intensity_grid(i,j,:); Grid_value=Grid_value(:); Intensity_value(i,j)=sum(Grid_value);   
        Grid_value2=Intensity_counts(i,j,:); Grid_value2=Grid_value2(:); Counts_value(i,j)=sum(Grid_value2);
    end
end
Intensity_value_2=Intensity_value;Intensity_value_2=rot90(Intensity_value_2,1);Intensity_value_2=flipud(Intensity_value_2);
Intensity_counts_2=Counts_value;Intensity_counts_2=rot90(Intensity_counts_2,1);Intensity_counts_2=flipud(Intensity_counts_2);

% save the output matrix


%% Step 2: calculate the CHRI ....

% Decomposite of the CHRI 1990-2020 & 2021-2050
clear

% load country and age group
population_65=readtable(' \World population 1990-2050 by age group\population age group 1990-2020.xlsx');
country_code=population_65.ISO3Alpha_code; country_code=string(country_code); pop_year=population_65.Year; pop_prop_total=population_65.Total; 
pop_prop_1=population_65.x50_69_1; pop_prop_2=population_65.x70__2;  pop_prop_3=pop_prop_1+pop_prop_2; 
country_id=readtable(' \Basic input files\Country ID.xlsx');
isocode=country_id.ISOCODE; isocode=string(isocode); country_value=country_id.Value_;
country_identifier=imread(' \Basic input files\national_identifier_grid_0.25deg1.tif');
country_identifier(country_identifier>1000)=0; Matrix_1=zeros(1,1440); country_identifier=[country_identifier; Matrix_1];

year_list=1990:2020;  

for ii=1:31

% Hazard,  normalized ETAI
filepath=' \CWs and HWs\temp_extremes_year*.mat';
filepath=dir(filepath); filepath2=[filepath(ii).folder '\' filepath(ii).name]; data=load(filepath2); data=data.data;
data=imresize(data,[720 1440]); input=data; input(isnan(input))=0;  input(input<=0)=NaN; input_1=input;  
top_thresh=prctile(input_1(:), 99); bot_thresh=prctile(input_1(:), 1); input_1(input_1>top_thresh)=top_thresh; input_1(input_1<bot_thresh)=bot_thresh;
haz_matrix= (input_1-min(min(input_1)))/(max(max(input_1))-min(min(input_1))); haz_matrix(haz_matrix<0.1)=NaN;

% Vulnerability, normalized amount of aging population     
age_population=load([' \Age structure 50+ 70+\CHRI_over50_31yrs\total_age_50_above_', num2str(year_list(ii)), '.mat']);  age_matrix_2=age_population.data; 
C=2; A=1;
vul_matrix = (log(A*age_matrix_2 + C) - log(min(min(A*age_matrix_2)) + C)) / (log(max(max(A*age_matrix_2)) + C) - log(min(min(A*age_matrix_2)) + C));  % or select: ln
vul_matrix(vul_matrix<=0.2)=NaN; 

% Exposure, normalized total population density  721*1440                       
filepath=' \Grided population density 1990-2020\pop_density_*.mat';
filepath=dir(filepath); filepath2=[filepath(ii).folder '\' filepath(ii).name]; data=load(filepath2); data=data.data; Matrix_1=zeros(1,1440); 
popul_den=[data; Matrix_1];  exp_matrix=(popul_den-min(min(popul_den)))/(max(max(popul_den))-min(min(popul_den))); exp_matrix(exp_matrix==0)=NaN; 

% determination and selection of w1 w2 w3 w4; refer to Fig. S9 in the submitted manuscript
EHRI_1 = w1*haz_matrix+w2*vul_matrix+w3*ada_matrix+w4*exp_matrix; 

% save files
end


%%  Step 3: tests ....

% trend and slope of time series ATEI and CHRI; grid by grid, global
input_data=total_ehri;
trend_result=zeros(720, 1440);   slope_result=zeros(720, 1440);
for i=1:720
    for j=1:1440 
        temp_vector=[]; temp_vector=input_data(i,j,:);  temp_vector= temp_vector(:); % temp_vector=temp_vector(temp_vector>0);  % Mann-Kendall test
        temp_vector2=temp_vector(~isnan(temp_vector));  kk1=isempty(temp_vector2); kk1=double(kk1); year2=1:length(temp_vector2); year2=year2';
        if kk1==0  && length(temp_vector2) >10          
            [H,p_value]=Mann_Kendall(temp_vector2,0.05);  %  [H,p_value]=Mann_Kendall(V,alpha); H=1 reject, H=0 fail to reject
            trend_result(i,j) = p_value; 
            slope_result(i,j) = Theil_Sen_Regress(year2,temp_vector2); % slope trend
        else
            trend_result(i,j) = NaN;  slope_result(i,j) = NaN;  
        end

    end
end
 result_1=slope_result; n1=numel(find(~isnan(result_1)));
 result_1(result_1<=0)=NaN; n2=numel(find(~isnan(result_1))); 

% Inequality in major countries-Gini

country_temp_order=sortrows(GDP_Temp2, 2); country_gdp_order=sortrows(GDP_Temp2, 1); gdp_temp_order_group_2=[];
for kk=1: size(GDP_Temp2,1)
    GDP_Temp_group=GDP_Temp2(kk,:);
    temp_order_row=find(country_temp_order(:,1)==GDP_Temp_group(1) & country_temp_order(:,2)==GDP_Temp_group(2)); temp_order_row=mean(temp_order_row);
    gdp_order_row=find(country_gdp_order(:,1)==GDP_Temp_group(1) & country_gdp_order(:,2)==GDP_Temp_group(2)); gdp_order_row=mean(gdp_order_row);
    gdp_temp_order_group_2(kk,:)=[temp_order_row gdp_order_row];
end
GDP_order_of_temp2=[gdp_temp_order_group_2(:,2) gdp_temp_order_group_2(:,1)]; GDP_orderdif_22=sortrows(GDP_order_of_temp2, 1); GDP_orderdif_32=GDP_orderdif_22(:,2);
x3_axis=[]; y3_axis=[]; 
for i=1:size(GDP_orderdif_32,1)
    x3_axis(i)=i/(size(GDP_orderdif_32,1)); 
    y3_axis(i)=sum(GDP_orderdif_32(1:i))/sum(GDP_orderdif_32);  
end
INEQ=(sum(y3_axis)-sum(y2_axis))/sum(y2_axis); 

%  slope of time series ETAI and EHRI; area by area; average value

clear

% function[H,p_value]=Mann_Kendall(V,alpha)
% V=reshape(V,length(V),1); 
% alpha = alpha/2; %
% n=length(V); 
% i=0; j=0; S=0; 
% for i=1:n-1
%    for j= i+1:n 
%       S= S + sign(V(j)-V(i)); 
%    end
% end
% VarS=(n*(n-1)*(2*n+5))/18;
% StdS=sqrt(VarS); 
%  %%% Note: ties are not considered 
% if S >= 0
%    Z=((S-1)/StdS)*(S~=0);
% else
%    Z=(S+1)/StdS;
% end
% p_value=2*(1-normcdf(abs(Z),0,1)); %% Two-tailed test 
% pz=norminv(1-alpha,0,1); 
% H=abs(Z)>pz; 
% return 

input=temp_extremes_aging_level;  [H,p_value]=Mann_Kendall(input(:,1),0.05);   trend_result_1 = p_value;
