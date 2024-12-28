function [out] = CENTERLINE(input)

    faces = input.faces;
    nodes = input.nodes;
    id_wall = input.id_wall;
    Laplace = input.Laplace;
        
    % Waitbar
    idum = 0;    
    estep = ceil(length(id_wall)*0.1);
    h = waitbar(0,{['Calculating centerline ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});
    
    centerline = zeros(length(id_wall),5);
    for n=1:length(id_wall)

        [cutpos,~,~,~] = qmeshcut(faces,nodes,Laplace,Laplace(id_wall(n)));
        cutpos = unique(cutpos,'rows');

        centerline(n,:) = [mean(cutpos), Laplace(id_wall(n)),n];
        
        idum = idum + 1;
        if mod(idum, estep) == 0
            waitbar(idum/ length(id_wall),h,{['Calculating centerline ... '],[num2str(idum),' nodes have been analyzed of a total ',num2str(length(id_wall)),'.']});
        end
    end

    centerl_out = sortrows(centerline,4);
    
    S2 = ceil((length(centerl_out(:,1))*0.05)); 
    idx = mod(S2,2)<1; 
    S2(idx) = S2(idx)+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if S2<3
        S2 = 3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out.centerline = centerl_out(:,1:3);
    out.centerline(:,1) = sgolayfilt(centerl_out(:,1),1,S2);
    out.centerline(:,2) = sgolayfilt(centerl_out(:,2),1,S2);
    out.centerline(:,3) = sgolayfilt(centerl_out(:,3),1,S2);

    out.centerline_lapid = centerl_out(:,4);
    out.centerline_n_id_wall = centerl_out(:,5);
 
    out.h=h;
    
    