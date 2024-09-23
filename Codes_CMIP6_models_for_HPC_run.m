clear
% model example for: ACCESS-CM2 
filepath1=' /CMIP6 2m air T data/SSP2-4.5 daily T max ACCESS-CM2 20210101-20501231.nc'; 
data_temp=ncread(filepath1,'tasmax');  data_temp=single(data_temp);  
indices=[1155	2616	4077	5538	6999	8460	9921];  data_temp(:,:,indices)=[]; 
  for ii=1:size(data_temp,3)
      Tmax=imresize(data_temp(:,:,ii),[1440 721],'bilinear');
      Tmax_Everyday_2(:,:,ii)=single(Tmax);   
  end
clearvars -except Tmax_Everyday_2  

 winter_threshold_1=load(' Daily mean T Threshold 9 models/ACCESS-CM2 winter_threshold.mat'); winter_threshold=winter_threshold_1.winter_threshold;
 summer_threshold_1=load(' /Daily mean T Threshold 9 models/ACCESS-CM2 summer_threshold.mat'); summer_threshold=summer_threshold_1.summer_threshold;
 winter_threshold=single(winter_threshold); summer_threshold=single(summer_threshold);
 data=load(' /Land_Water_01_Global_Grid_0.25deg.mat'); land_water=data.data; Matrix5_1=zeros(1,1440); 
 land_water=[Matrix5_1; land_water]; land_water=[land_water(:,721:1440) land_water(:,1:720)]; land_water=rot90(land_water,3);  
 
 
 Year_num=29;  start_day=date1; duration=151;   %  date2
 HW_central=zeros(1440, 721,duration, Year_num, 'single'); HW_neighbor=zeros(1440, 721,duration, Year_num, 'single'); % heat waves
 CW_central=zeros(1440, 721,duration, Year_num, 'single'); CW_neighbor=zeros(1440, 721,duration, Year_num, 'single'); % cold spells
 for i=2:1439
     for j=2:720  
        day_count=0;   
        if land_water(i,j-1)==1   
            for d=start_day:start_day+duration-1  %calendar day May-August 121-243; November-February
                cout_num=0; day_count=day_count+1; Sliding_cell_total=[];
        
                for m=1:Year_num   % total 30 years             
                  for n=d-3+365*(m-1):d+3+365*(m-1)  % total 7 days  sliding window     
                     cout_num=cout_num+1;
                     Neighbor_grid_cells=Tmax_Everyday_2(i-1:i+1,j-1:j+1,n); 
                     Sliding_cell_total(:,cout_num)=Neighbor_grid_cells(:);  
                  end        
                end    
              Threshold_h_90th=prctile(Sliding_cell_total(:),90); Threshold_c_10th=prctile(Sliding_cell_total(:),10); 
                 for mm=1:Year_num   % total 30 years             
                    Central_grid=Tmax_Everyday_2(i,j,d+365*(mm-1));                
                    Neighbor_grid_cells=[];  Neighbor_grid_cells=Tmax_Everyday_2(i-1:i+1,j-1:j+1,d+365*(mm-1));
                    threshold_sum=summer_threshold(i,j,mm);  threshold_win=winter_threshold(i,j,mm);
                    if Central_grid > threshold_sum                    
                        Intensity=abs(Central_grid/Threshold_h_90th); HWs_count=length(Neighbor_grid_cells>Threshold_h_90th); %  heat waves
                        HW_central(i, j, day_count, mm)=single(Intensity); HW_neighbor(i, j, day_count, mm)=single(HWs_count);   
                    end
                    if Central_grid < threshold_win     
                        Intensity=abs(Threshold_c_10th/Central_grid); CWs_count=length(Neighbor_grid_cells<Threshold_c_10th); %  cold spells
                        CW_central(i, j, day_count, mm)=single(Intensity); CW_neighbor(i, j, day_count, mm)=single(CWs_count);
                    end
                end   
         
            end
        else
            continue % loop to next grid
        end
    end
end

Intensity_grid=zeros(1439, 720,Year_num, 'single'); Intensity_counts=zeros(1439, 720,Year_num, 'single');
for i=2:1439
    for j=2:720
        for k=1:Year_num 
            Single_year_series= HW_central(i,j,:,k);   Single_year_neibour= HW_neighbor(i,j,:,k);
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
filename_1=' /Daily max T 9 models/ACCESS-CM2 intensity 2021-2050 winter HWs.mat'; save(filename_1, 'Intensity_grid');
filename_2='/ /Daily max T 9 models/ACCESS-CM2 counts 2021-2050 winter HWs.mat'; save(filename_2, 'Intensity_counts');

Intensity_grid=zeros(1439, 720,Year_num, 'single'); Intensity_counts=zeros(1439, 720,Year_num, 'single');
for i=2:1439
    for j=2:720
        for k=1:Year_num 
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
filename_1=' /Daily max T 9 models/ACCESS-CM2 intensity 2021-2050 winter CWs.mat'; save(filename_1, 'Intensity_grid');
filename_2=' /Daily max T 9 models/ACCESS-CM2 counts 2021-2050 winter CWs.mat'; save(filename_2, 'Intensity_counts');
