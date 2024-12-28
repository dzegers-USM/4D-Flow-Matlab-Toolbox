function [out] = LENGTH_AND_CURVATURE(input)

    centerline = input.centerline;
    n_id_wall = input.centerline_n_id_wall;
    id_wall = input.id_wall;
    Laplace = input.Laplace;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vessel distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dis_bet_p_cent = sqrt(sum((centerline(1:end-1,:)-centerline(2:end,:)).^2,2));
    d = [0;cumsum(dis_bet_p_cent)]; % length of the centerline for each point
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change point of the centerline to 2.5 mm equal spacing %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = waitbar(0,'Wait adjusting centerline for curvature quantification ... ');
        
    if d(end)>30
        
        num_p = round(d(end)/2.5);
        cc = centerline([1:end 1],:)'/1000;
        ff = fnplt(cscvn(cc(:,1:end-1)))';
        ff = unique(ff,'rows','stable');
        tt = linspace(0,1,num_p);
        pt = interparc(tt,ff(:,1),ff(:,2),ff(:,3),'pchip');

        centerline_new_1 = pt; % centerline equal spacing to 2.5 mm in meters
        gg = sqrt(sum((centerline_new_1(1:end-1,:)-centerline_new_1(2:end,:)).^2,2)); % spacing between each point of the new centerline
        tt = double([0;cumsum(gg)]); % total length of the centerline in m

        x_1d = gradient(centerline_new_1(:,1),gg(1));
        y_1d = gradient(centerline_new_1(:,2),gg(1));
        z_1d = gradient(centerline_new_1(:,3),gg(1));

        x_2d = gradient(x_1d,gg(1));
        y_2d = gradient(y_1d,gg(1));
        z_2d = gradient(z_1d,gg(1));

        % curvature to the new centerline and interpolation
        curvature_cent = sqrt((z_2d.*y_1d - y_2d.*z_1d).^2 + (x_2d.*z_1d - z_2d.*x_1d).^2 + (y_2d.*x_1d - x_2d.*y_1d).^2)./((x_1d.^2 + y_1d.^2 + z_1d.^2).^(3/2)); 
        cc = feval(fit(tt,curvature_cent,'SmoothingSpline','SmoothingParam',0.999997),tt);
%         cc_n = interp1(tt*1000,cc,linspace(tt(1),tt(end),size(centerline,1)),'pchip');
        cc_n = interp1(tt*1000,cc,d,'pchip');

        % Waitbar
        idum = 0;    
        estep = ceil(length(id_wall)*0.1);
        close(h)
        h = waitbar(0,{['Calculating length of the vessel and curvature ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});

        length_vessel = zeros(size(Laplace));
        curvature = zeros(size(Laplace));
        for n = 1:length(id_wall)

            [r,~,~] = find(n_id_wall==n);
            length_vessel(id_wall(n),1) = d(r);
            curvature(id_wall(n),1) = cc_n(r);

            idum = idum + 1;
            if mod(idum, estep) == 0
                waitbar(idum/ length(id_wall),h,{['Calculating length of the vessel and curvature ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});
            end
        end

        out.length_vessel = length_vessel;
        out.curvature = curvature;
        out.reduced_centerline = centerline_new_1*1000; % reduced centerline in mm;
        out.h = h;
        
    else
        
        waitfor(warndlg('To calculate the curvature the vessel need to be larger than 3 cm!','Warning'))
        
        % Waitbar
        idum = 0;    
        estep = ceil(length(id_wall)*0.1);
        close(h)
        h = waitbar(0,{['Calculating length of the vessel and curvature ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});

        length_vessel = zeros(size(Laplace));
        for n = 1:length(id_wall)

            [r,~,~] = find(n_id_wall==n);
            length_vessel(id_wall(n),1) = d(r);

            idum = idum + 1;
            if mod(idum, estep) == 0
                waitbar(idum/ length(id_wall),h,{['Calculating length of the vessel and curvature ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});
            end
        end

        out.length_vessel = length_vessel;
        out.curvature = zeros(size(Laplace));
        out.reduced_centerline = centerline; % reduced centerline in mm;
        out.h = h;
        
    end
    
    

