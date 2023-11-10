function [nt,xstream2,ystream2,xs,ys,vvs,vx_srf,vy_srf,xi_all,yi_all,Xs,Ys,Fvxs,Fvys] = DDTool_ISMIP6(dist,...
    np,startpoints_type,X,Y,u,v,thk,ix,iy,clim_in, direction,backS,varargin)
% This function generates streamlines from which catchments can be delineated. It takes as input: 
% dist = distance between points in the streamline, 
% np = maximun number of points in a streamline
% startpoint_types: string input {'clim'|'specific_coord'| 'specific_index'| 'sensitivity'},
% X,Y: grid coordinate vectors from the given ISMIP6 ice flow model
% u,v,thk: modeled easting and northing velocity fields and ice thickness from a given ice flow model.
% ix,iy: streamline starting point, can eaither be a coordinate, an index, or empty [] depending on the
% chosen 'startpoint_types'. 
% clim_in: minimum velocity of gridpoint from which a streamline can be started, (optional, 
% only required if 'startpoint_types'= 'clim')
% direction: string describing if flowline show be calculated along or
% against flow direction, or in both directions i.e. {'bwd'|'fwd'|'fwd_bwd'}
% backS = optional string 'backS', if present, function calculates back streamline. The Function 
% is called seperately for this.

% The function outputs: 
% nt: number of time slices
% xstream2, ystream2: coordinates of all the started streamlines,
% xs, ys; coordinates of the grid from the smaller working region, 
% vvs,vx_srf,vy_srf
% xi_all, yi_all: coordinates at which flow lines are started if 'startpoint_types'='clim' 
% Xs,Ys: Grid vectors reduced to working region
% Fvxs,Fvys: Interpolated surfaces of the easting and northing velocity components

% Note, this tool was made specifically to work with the ISMIP6 model output files.

%% 1. Change the potential empty first layer 

if sum(~isnan(u(:,:,1)),'all') == 0
    disp('no velocity values')
    u(:,:,1) = u(:,:,2); 
    v(:,:,1) = v(:,:,2); 
end

%% 2. Decreasing the working region 

% define the number of timeslices 
dvel = size(u); % dimension of the velocity field
nt = dvel(end); % number of time slices

% define the working region encompassing the entire JI catchement 
ixs = find(X>-250000 & X<400000);
iys = find(Y<-2*10^6 & Y>-2.4*10^6);

Xs = X(ixs); 
Ys = Y(iys); 
[ys,xs] = meshgrid(Ys,Xs);

% define new matrix dimensions 
lx = length(Xs); 
ly = length(Ys);

% limit the input maps to the working region 
vx_srf = u(ixs,iys,:);
vy_srf = v(ixs,iys,:);
thk = thk(ixs,iys,:);

% calculate the total velocity magnitude from the components
vvs = sqrt(vx_srf.^2+vy_srf.^2);

%% 3. Calculating flow direction, correcting to be w respect to north i.e. azimuth

theta_r = zeros(lx,ly);   

for k = 1:nt
for i = 1:lx
    for j= 1:ly
        if vx_srf(i,j,k)>0                          
            if vy_srf(i,j,k)>0  
                theta_r(i,j,k) = atand(vx_srf(i,j,k)/vy_srf(i,j,k));      % I
            else            
                theta_r(i,j,k) = atand(vy_srf(i,j,k)/vx_srf(i,j,k))+90;   % IV
            end
        else                              
            if vy_srf(i,j,k)>0  
                theta_r(i,j,k) = atand(vy_srf(i,j,k)/vx_srf(i,j,k))+270;  % II
            else                                
                theta_r(i,j,k) = atand(vx_srf(i,j,k)/vy_srf(i,j,k))+180;  % III     
            end
        end
    end
end
end

%% 4. Calculating the streamlines 

for k = 1:nt
%% 4.1 Creating interpolated surfaces 
    Fvxs = griddedInterpolant(xs,ys,vx_srf(:,:,k)); %vx_srf(:,:,k)
    Fvys = griddedInterpolant(xs,ys,vy_srf(:,:,k));
    Fthetar = griddedInterpolant(xs,ys,theta_r(:,:,k));
    Fvvs = griddedInterpolant(xs,ys,vvs(:,:,k)); % test field

%% 4.2 Back Streamline  

if strcmp(direction,'bwd') && nargin > 12 
    disp(nargin)
    % start point of back stream line is the end point of outermost streamline 
    sprintf('now calculate the back end stream line, closing the domain')
    
    % defining the starting point of the back streamline 
    xi = ix;
    yi = iy(1);
    bxstream2(1,1)=ix; 
    bystream2(1,1)=iy(1); 

    for i = 1:np % go until I hit the second line!
        
        % calculate the timestep size required to move the specified
        % distance 'dist' with the velocity at the quiry point.
        btime(:,i) = dist/Fvvs((bxstream2(:,i)),(bystream2(:,i)));

        % calculating the first streamline point coordinates
        bdxstream2(:,i) = -1*Fvxs((bxstream2(:,i)),(bystream2(:,i)))*btime(:,i); 
        bdystream2(:,i) = -1*Fvys((bxstream2(:,i)),(bystream2(:,i)))*btime(:,i);

        % update the streamline coordinate 
        bxstream2(:,i+1) = (bxstream2(:,i)) + bdxstream2(:,i);
        bystream2(:,i+1) = (bystream2(:,i)) + bdystream2(:,i); %error near front
        
        % stop streamline calculation if it extents beyond the other
        % boundary.
        if bystream2(:,i+1) > iy(2)
            fprintf('y is larger than end y point')
            break;
        end

        % retrieving the flow direction values at the updated points
        btheta_loop2(:,i+1) = Fthetar(bxstream2(:,i+1),bystream2(:,i+1));
    end

    % defining output variables
    xstream2 = bxstream2; ystream2 = bystream2;

    return
end

%% 4.3 Normal streamline calculations   

switch startpoints_type

    case 'clim' % find streamline starting gridpoints based on velocity limit
        clim = clim_in; % setting the velocity limit were to start gridpoints
        
        % defining search region near the front of the ice sheet.  
        [it1,it2] = find(ys<-2.24*10^6 & ys>-2.295*10^6 & xs>-2.23*10^5 & -1.5*10^5>xs);
        
        kk = find(vvs(unique(it1),unique(it2),k)>clim); % no of gridpoints to start streamlines from 
        
        % re-sizing the coordinates to only be within the small search
        % region near the front, otherwise the kk index wont match coordinate indecies  
        xs2 = xs(unique(it1),unique(it2)); 
        ys2 = ys(unique(it1),unique(it2));
        ys3 = ys2(:); % reshaping to cloumn vector to match the index kk vector
        xs3 = xs2(:);
        yi = ys3(kk);
        xi = xs3(kk);
        
        yi_all{k} = yi;
        xi_all{k} = xi;

    case 'specific_coord' % Input to function manually defining the streamline starting points
        xi = ix;
        yi = iy;
        
        % need to output these in the same format as the clim, though
        % should be the same for all k layers!
        yi_all{k} = yi;
        xi_all{k} = xi;
        
    case 'specific_index' % selecting already existing gridpoints to start flowlines from 

        % starting point
        xi = Xs(ix)';
        yi = Ys(iy)';

    case 'sensitivity' % starting several streamlines at various y- coordinates at given x-coordinate
        if length(ix)>1 
            fprintf('only one x-index value required')
            return
        elseif length(iy)>4
            fprintf('only 4 y-indecies required')
            return
        end
        l1 = (iy(2)-iy(1))*3;
        l2 = (iy(4)-iy(3))*3;
        yi = [linspace(Ys(iy(1)),Ys(iy(2)),l1) linspace(Ys(iy(3)),Ys(iy(4)),l2)];
        xi = ones(1,length(yi)).*Xs(ix);
end


switch direction
    case {'fwd','bwd'} % The first part is the same for both cases therefore collect it 
        for j = 1:length(yi)
            
            % setting up the cell and defining the first entrance for each stream line as the initial point value
            xstream2{k,j}=xi(j);
            ystream2{k,j}=yi(j);

            switch startpoints_type
                case {'specific_coord','clim'}
                    theta_loop2{j} = Fthetar(xstream2{j},ystream2{j}); 
                case 'specific_index'
                    thetai{j}=theta_r(ix(j),iy(j));
                    theta_loop2{j}=thetai{j};
                case 'sensitivity'
                    theta_loop2{j} = Fthetar(xstream2{j},ystream2{j});
            end

            for i = 1:np-1 
                % stop streamlines exiting the smaller working region
                if xstream2{k,j}(i)>xs(end)-dist || ystream2{k,j}(i)>ys(end)-dist
                    fprintf('Coordinate out of bounds \n')
                    break;
                end

                vel_min = 0.5; % [m/yr]
                % stop streamline when reaching specific ice flow velocity
                if Fvvs(xstream2{k,j}(i),ystream2{k,j}(i))< vel_min 
                    break
                end

                switch direction
                    case 'bwd'
                        % reversing the flow direction 180 degrees to 'walk' backwards along ice path
                        if (theta_loop2{j}(i)) < 180
                            theta_loop2{j}(i) = (theta_loop2{j}(i))+180; 
                        else
                            theta_loop2{j}(i) = (theta_loop2{j}(i))-180;
                        end

                        time{k,j}(i) = dist/Fvvs((xstream2{k,j}(i)),(ystream2{k,j}(i)));

                        % calculating the first streamline point coordinates
                        dxstream2{k,j}(i) = -1*Fvxs((xstream2{k,j}(i)),(ystream2{k,j}(i)))*time{k,j}(i);
                        dystream2{k,j}(i) = -1*Fvys((xstream2{k,j}(i)),(ystream2{k,j}(i)))*time{k,j}(i);

                        % Updating the streamline
                        xstream2{k,j}(i+1) = (xstream2{k,j}(i)) + dxstream2{k,j}(i);
                        ystream2{k,j}(i+1) = (ystream2{k,j}(i)) + dystream2{k,j}(i);

                        % retrieving the flow direction values at the updated points
                        theta_loop2{j}(i+1) = Fthetar(xstream2{k,j}(i+1),ystream2{k,j}(i+1));

                    case 'fwd'
                        
                        time{k,1,j}(i) = dist/Fvvs((xstream2{k,1,j}(i)),(ystream2{k,1,j}(i)));

                        % calculating the first streamline point coordinates
                        dxstream2{k,1,j}(i) = Fvxs((xstream2{k,1,j}(i)),(ystream2{k,1,j}(i)))*time{k,1,j}(i); 
                        dystream2{k,1,j}(i) = Fvys((xstream2{k,1,j}(i)),(ystream2{k,1,j}(i)))*time{k,1,j}(i);

                        xstream2{k,1,j}(i+1) = (xstream2{k,1,j}(i)) + dxstream2{k,1,j}(i);
                        ystream2{k,1,j}(i+1) = (ystream2{k,1,j}(i)) + dystream2{k,1,j}(i);%error near front

                end
            end
        end

    case 'fwd_bwd' % streamline in both direction from starting point    
        for j = 1:length(yi)    
        
        % need to save 2 streamlines for each startingpoints.
        xs_fwdbwd{k,j}(1:2)=xi(j);
        ys_fwdbwd{k,j}(1:2)=yi(j);

        switch startpoints_type
            case {'specific_coord','clim'}
                theta_loop2{k,j} = Fthetar(xs_fwdbwd{k,j},ys_fwdbwd{k,j}); 
            case 'specific'
                thetai{j}=theta_r(ix(j),iy(j));
                theta_loop2{j}=thetai{j};
            case 'sensitivity'
                theta_loop2{k,j} = Fthetar(xs_fwdbwd{k,j}(2),ys_fwdbwd{k,j}(2));
        end             

        for i = 1:np-1

            % stop streamlines exiting the smaller working region
            if xs_fwdbwd{k,j}(i,2)>xs(end)-dist || ys_fwdbwd{k,j}(i,2)>ys(end)-dist
                fprintf('Coordinate out of bounds \n')
                break;
            end

            % stop streamline when reaching specific ice flow velocity
            if Fvvs(xs_fwdbwd{k,j}(i,2),ys_fwdbwd{k,j}(i,2))< 0.2
                break
            end

            % time calc
            time_fwd{k,j}(i) = dist/Fvvs((xs_fwdbwd{k,j}(i,1)),(ys_fwdbwd{k,j}(i,1)));
            time_bwd{k,j}(i) = dist/Fvvs((xs_fwdbwd{k,j}(i,2)),(ys_fwdbwd{k,j}(i,2)));

            % backwards angle calc
            if (theta_loop2{k,j}(i)) < 180 
                theta_loop2{k,j}(i) = (theta_loop2{k,j}(i))+180; 
            else
                theta_loop2{k,j}(i) = (theta_loop2{k,j}(i))-180;
            end

            % calculating the first streamline point coordinates
            dx_fwd{k,j}(i) = Fvxs((xs_fwdbwd{k,j}(i,1)),(ys_fwdbwd{k,j}(i,1)))*time_fwd{k,j}(i);
            dy_fwd{k,j}(i) = Fvys((xs_fwdbwd{k,j}(i,1)),(ys_fwdbwd{k,j}(i,1)))*time_fwd{k,j}(i);

            dx_bwd{k,j}(i) = -1*Fvxs((xs_fwdbwd{k,j}(i,2)),(ys_fwdbwd{k,j}(i,2)))*time_bwd{k,j}(i);
            dy_bwd{k,j}(i) = -1*Fvys((xs_fwdbwd{k,j}(i,2)),(ys_fwdbwd{k,j}(i,2)))*time_bwd{k,j}(i);

            xs_fwdbwd{k,j}(i+1,1) = (xs_fwdbwd{k,j}(i,1)) + dx_fwd{k,j}(i);
            ys_fwdbwd{k,j}(i+1,1) = (ys_fwdbwd{k,j}(i,1)) + dy_fwd{k,j}(i);

            xs_fwdbwd{k,j}(i+1,2) = (xs_fwdbwd{k,j}(i,2)) + dx_bwd{k,j}(i);
            ys_fwdbwd{k,j}(i+1,2) = (ys_fwdbwd{k,j}(i,2)) + dy_bwd{k,j}(i);

            % retrieving the flow direction values at the updated points
            theta_loop2{k,j}(i+1) = Fthetar(xs_fwdbwd{k,j}(i+1,2),ys_fwdbwd{k,j}(i+1,2));
        end
        end 

        % save to the common output variable name
        xstream2 = xs_fwdbwd;
        ystream2 = ys_fwdbwd;
end

clear Fthetar Fvvs
end

end 