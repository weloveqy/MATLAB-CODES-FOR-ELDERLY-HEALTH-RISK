
%% per hour ERA5 LST to daily max min mean Temperature
clear

% ERA5 air temperature 0.25 degree
listing = dir(' \ERA5 2m air temperature\t2m_*24h.nc');
list_size=size(listing,1);
for l=1:list_size
    name=listing(l).name; filepath=[' \ERA5 2m air temperature\' name];
    Info=ncinfo(filepath); Size=Info.Dimensions(3).Length; Size=Size/24; Year=filepath(44:47);
    Tmax_matrix=zeros(1440,721); Tmin_matrix=zeros(1440,721); Tmean_matrix=zeros(1440,721);

    for k=1:Size
     data_t2m=ncread(filepath,'t2m',[1,1,24*(k-1)+1],[1440,721,24]);
        for i=1:1440
         for j=1:721
         data_t2m_24h=data_t2m(i,j,1:24);   data_t2m_24h=squeeze(data_t2m_24h(1,1,:)); plot(data_t2m_24h)
         data_t2m_24h=squeeze(data_t2m_24h);
         Tmax=max(data_t2m_24h);Tmin=min(data_t2m_24h);Tmean=mean(data_t2m_24h);
         Tmax_matrix(i,j)=Tmax; Tmin_matrix(i,j)=Tmin; Tmean_matrix(i,j)=Tmean;
         end
        end
    filename_1=[' \ERA5 2m air temperature\Year_' num2str(Year) '_Day_' num2str(k) '_Tmax_Tmin_Tmean.mat'];
    save(filename_1, 'Tmax_matrix', 'Tmin_matrix', 'Tmean_matrix');
    end

end
