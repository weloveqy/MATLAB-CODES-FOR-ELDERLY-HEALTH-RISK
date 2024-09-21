
%% annual mean temperature, 75th percentile summer, 25th percentile winter

Year_num=33; summer_threshold=zeros(1439,721,Year_num); winter_threshold=zeros(1439,721,Year_num);
cycle_count=0;
data=load(' \Land_Water_01_Global_Grid_0.25deg.mat');
land_water=data.data;  
land_water=flipud(land_water); land_water=rot90(land_water,3);  

for i=2:1439
    for j=2:720
       cycle_count=cycle_count+1; cycle_count
        if land_water(i,j-1)==1           
            for m=1:Year_num   % total 33 years                  
                slide_window=Tmean_Everyday(i-1:i+1,j-1:j+1,365*(m-1)+1:365*m); 
                sliding_cells=slide_window(:);         
                Q1=prctile(sliding_cells,25); Q3=prctile(sliding_cells,75); 
                summer_threshold(i,j,m)=Q3; winter_threshold(i,j,m)=Q1;
            end
        else
            continue
        end
    end
end