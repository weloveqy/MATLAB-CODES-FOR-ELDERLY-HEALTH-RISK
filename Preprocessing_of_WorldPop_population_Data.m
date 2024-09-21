%% Preprocessing of WorldPop population Data 

%% read and transfer the resolution 
clear

filepath=' \global_*_2020_1km.tif'; filepath=dir(filepath); 

for ii=1:length(filepath)
    num=ii; num
    filepath2=[filepath(ii).folder '\' filepath(ii).name]; data=imread(filepath2);
    population_data=data; population_data(population_data<0)=0; age_matrix=zeros(624,1440);   
    for i=1:30:18720
        for j=1:30:43200
            m=(i-1)/30+1; n=(j-1)/30+1;
            matrix=population_data(i:i+29,j:j+29);  sum_1=sum(matrix(:)); age_matrix(m,n)=sum_1;                   
        end
    end
    Matrix_1=zeros(24,1440); Matrix_2=zeros(72,1440); age_matrix_2=[Matrix_1; age_matrix; Matrix_2]; age_matrix_2(age_matrix_2<0)=0;
    age_50_above(:,:,ii)=age_matrix_2;
end

for i=1:720
    for j=1:1440
        matrix= age_50_above(i,j,:);  matrix=matrix(~isnan(matrix));
        total_age_50_above(i,j)=sum(matrix(:));
    end
end
data=total_age_50_above;
filename_1=' \total_age_50_above.mat'; save(filename_1, 'data'); 

%% produce age structures 1990-1999 2020-2050 % missing data; deal with data gap
clear

data=load(' \Age structure 50+ 70+\CHRI_over50_31yrs\total_age_50_above_2020'); data1=data.data; % imagesc(data1)

% load country and age group
population_65=readtable(' \World population 1990-2050 by age group\population age group 2021-2050.xlsx');
country_code=population_65.ISO3Alpha_code; country_code=string(country_code); pop_year=population_65.Year; 
pop_prop_50=population_65.x50_69+population_65.x70__1;  pop_prop_70=population_65.x70__1;

country_id=readtable(' \Basic input files\Country ID.xlsx');
isocode=country_id.ISOCODE; isocode=string(isocode); country_value=country_id.Value_;
country_identifier=imread(' \Basic input files\national_identifier_grid_0.25deg1.tif');
country_identifier(country_identifier>1000)=0;  % imagesc(country_identifier);

year_list=2021:2050;  % year_list=1990:1999;  
for ii=1:30
aged_group_grid=zeros(720,1440);
for i=1:720
    for j=1:1440
        value_1=country_identifier(i,j); 
        if value_1>0 
            row=find(country_value==value_1); value_2=isocode(row); row_2=find(country_code==value_2 & pop_year==year_list(ii));
            row_3=find(country_code==value_2 & pop_year==2020);
            if row_2>0              
                aged_prop=pop_prop_50(row_2,1); aged_prop_2=pop_prop_50(row_3,1);
                aged_group_grid(i,j)=data1(i,j)*aged_prop/aged_prop_2; 
            end
        else
            aged_group_grid(i,j)=NaN;  % imagesc(aged_group_grid);
        end
    end
end 

data=aged_group_grid; 
save([' \CHRI_over50_31yrs\total_age_50_above_', num2str(year_list(ii)), '.mat'], 'data');

end