% Delineating Catchments from ISMIP6 model outputs
clc; clear all; close all;

% Flowline function parameters 
np = 550;                       % number of points in the streamline
dist = 1000;                    % distance between points on the 'streamline'
clim = 10;                      % velocity limit starting streamline points from 
x_outerlines = -1.5*10^5;       % cut-off x-coordinate

%% Opening the netcdf files - updated 
% ISMIP6 save individual files for each variable, therefore make structure 
% containing the name of each file for each variable for download

addpath('Path\to\the\ISMIP6\Files')

% name strings to search for different comparisons:
EXPs = "_historical";

fp_ismip = 'C:\Users\njlk\Documents\ISMIP6_output_download'; % filepath for the ismip6 .nc files 
vxf = dir(fullfile(fp_ismip,'xvelsurf*.nc'));                % structure w filenames of vx
vyf = dir(fullfile(fp_ismip,'yvelsurf*.nc'));                % structure w filenames of vy
thkf = dir(fullfile(fp_ismip,'lithk*.nc'));                  % structure w filenames of thk

% Getting the common grid from PISM files since it contain x, not all files does!
% these variables should not change between the different files
X = ncread(vxf(end).name,'x'); % [m]
Y = ncread(vxf(end).name,'y'); % [m]

% Selecting only a subset of the files in the folder 
count = 1;
idx_name = [];
check = EXPs; 

for ii = 1:length(vxf)
    for jj = 1:length(check)
        if ~isempty(strfind(vxf(ii).name,check(jj)))
            idx_name(count) = ii;
            count = count+1;
        end
    end
end

vxf = vxf(idx_name);
vyf = vyf(idx_name);
thkf = thkf(idx_name);

%% delineating the catchment 

Legend=cell(length(vxf),1); % preallocating space for legend entries for later plotting
nfts_all = zeros(1,length(vxf)); 

for i = 1:length(vxf) % for each group participating in ISMPI6
close(figure(i));figure(i)
title(append(strrep(vxf(i).name(14:end-3),'_',' ')))

        disp(i)
        % retrieving the variables from a given group run
        vx = ncread(vxf(i).name,'xvelsurf').*(60*60*24*365);  % [m/s]--> [m/yr]
        vy = ncread(vyf(i).name,'yvelsurf').*(60*60*24*365);  % [m/s]--> [m/yr]
        thk = ncread(thkf(i).name,'lithk');                   % [m]
        
        % if only working with last time slice:
        vx = vx(:,:,end-1:end);
        vy = vy(:,:,end-1:end);
        thk = thk(:,:,end-1:end);
        
        % changing zero values to NaN
        vx(vx==0) = nan;
        vy(vy==0) = nan;
        
        % defining the run type
        type = 'clim';                  
        
        % calling the function calculating the stream lines
        [nt,bxs,bys,xs,ys,vvs,vx_s,vy_s,xi,yi,Xs,Ys,Fvxs,Fvys] = DDTool_ISMIP6(dist,np,...
            type,X,Y,vx,vy,thk,[],[],clim,'bwd');
                
        size_f = size(vx_s);
        nfts = size_f(3); % number of time slices in a given file with non-zero velocity maps
        nfts_all(i) = nfts;
        
        % saving all the shortened velocity fields for later plotting
        vvs_tot(:,:,i) = vvs(:,:,nfts);
        vxs_tot(:,:,i) = vx_s(:,:,nfts);
        vys_tot(:,:,i) = vy_s(:,:,nfts);
        
        for k = nfts  % Finding the catchment for only the final time slice i.e. 2015 I hope!
            
            % plotting the background velocity map
            figure(i); 
            pcolor(xs,ys,log(vvs(:,:,k))); hold on; shading flat;
            axis equal tight; caxis([-2 7]);
            title(['v_{clim} = ',num2str(clim)]);
            plot(x_outerlines+zeros(1,length(ys(1,:))),ys(1,:),':k');
            plot(xi{k},yi{k},'.k');     % plot gridpoints starting flowlines 
            
            maxy = [];
            miny = [];
            
            [numRows,numCols] = size(bxs);
            for j = 1:(numCols-sum(cellfun(@isempty,bxs(k,:))))
                
                % Plotting all the generated streamlines
                plot(bxs{k,j},bys{k,j},'-','Color', [.7 .7 .7]);
                
                % find the 2 outermost streamlines delineating the drainage domain:
                maxy(j) = max(bys{k,j}(bxs{k,j}<x_outerlines));
                miny(j) = min(bys{k,j}(bxs{k,j}<x_outerlines));
                
                % method 2 finding outermost streamlines, find and sum y-values
                % within a range
                xlimit = -1*10^5;
                ylimit = -2.27*10^6;
                
                uplines_y = bys{k,j}(bxs{k,j}<(-1*10^5) & bys{k,j}>-2.27*10^6);
                uplines_x = bxs{k,j}(bxs{k,j}<(-1*10^5) & bys{k,j}>-2.27*10^6);
                
                lUp(j) = length(uplines_y);
                
                y_range = bys{k,j}(bxs{k,j}<xlimit);
                y_range_max = bys{k,j}(bxs{k,j}<(-1*10^5) & bys{k,j}>(-2.27*10^6));
                ysum(j) = sum(y_range);
                
                l_yra(j) = length(y_range); % length/number of points within the range
                l_yra2(j) = length(y_range_max);
                
                % trying to interpolate to find outermost stream line
                if sum(isnan(bxs{k,j})) > 100 || length(bxs{k,j})<10
                    yq1(j) = nan;
                else
                    [ubxs,idxu,idxc]=unique(bxs{k,j});
                    ubys = bys{k,j}(idxu);
                    yq1(j) = interp1(ubxs,ubys,-1.6*10^5); 
                end
                
            end
            imin = max(maxy); disp(imin);
            imax = min(miny); disp(imax);
            
            Imax(k) = find(maxy==max(maxy));
            Imin(k) = find(miny==min(miny));
            
            % Plotting the found outermost streamlines from which the
            % catchment will be defined
            plot(bxs{k,Imax(k)},bys{k,Imax(k)},'*-w');
            plot(bxs{k,Imin(k)},bys{k,Imin(k)},'-*w');
            
            % Calling the function to calculate the back streamline, using the end
            % points of each of the outermost streamlines
            type = 'specific_coord'     % 'sensitivity'
            
            % back flowline starting from southern boundary
            iyb2(1) = bys{k,Imin(k)}(1,end); 
            ixb2 = bxs{k,Imin(k)}(1,end);    
            iyb2(2) = bys{k,Imax(k)}(1,end); 
                        
            % calling the function 
            if strcmp(vxf(i).name,'xvelsurf_GIS_JPL_ISSM_historical.nc') == 1
               
               % midpoint to northen boundary  
               iyb21(1) = -2.2225*10^6; 
               ixb21 = 2.84*10^5;          
               iyb21(2) = bys{k,Imax(k)}(1,end);  
                
               type = 'specific_coord'
               [nt21,bxs21,bys21] = DDTool_ISMIP6(dist,np,type,X,Y,vx,vy,thk,ixb21,...
                   iyb21,clim,'bwd','backS');
                
               % midpoint to southern boundary
                iyb22(1) = -2.242*10^6;  
                ixb22 = 2.86*10^5;     
                iyb22(2) = -2.242*10^6;  

                [nt22,bxs22,bys22] = sensitivity_DDtool_ismip6(dist,np,type,X,Y,vx,vy,thk,ixb22,...
                   iyb22,clim,'bwd','backS');
               
                bxs2 = [flip(bxs22(bys22>bys{k,Imin(k)}(1,end))) bxs21];
                bys2 = [flip(bys22(bys22>bys{k,Imin(k)}(1,end))) bys21]; 
                
            else
                [nt,bxs2,bys2] = DDTool_ISMIP6(dist,np,type,X,Y,vx,vy,thk,ixb2,iyb2,...
                    clim,'bwd','backS');                
            end
            
            % plotting the 'backstreamline' closing the catchment
            figure(i); hold on;
            plot(bxs2,bys2,'-m');
                
            ddx1{i,k} = [bxs{k,Imin(k)} bxs2(2:end) flip(bxs{k,Imax(k)})];
            ddy1{i,k} = [bys{k,Imin(k)} bys2(2:end) flip(bys{k,Imax(k)})];
                
        end
        
        figure(100+i);
        title(strrep(vxf(i).name(14:end-3),'_',' '))
        hold on
        plot(ddx1{i,2}./1000,ddy1{i,2}./1000,'Linewidth',1.2);
        plot(xi{k}./1000,yi{k}./1000,'.');
        axis equal tight
        xlabel('x_p [km]'); ylabel('y_p [km]')
        
        % calculate the area
        polyin = polyshape([ddx1{i,k}' ddy1{i,k}']); 
        A(i) = area(polyin);
            
end
