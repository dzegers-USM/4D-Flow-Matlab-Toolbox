function autoSelectSectionsLaplace(handles)
    handles.sec = str2num(get(handles.edit3,'String'));
    if handles.sec > 1
        % values of length that identify the section
        values_sec = linspace(0,max(handles.length_vessel),handles.sec+1);
        ID_nodes = [1:size(handles.nodes,1)]';
        
        % laplace node selected
        ID_selected = zeros(handles.sec + 1,1);
        for n = 2:2+(handles.sec-2)
            [~,I] = min(sqrt((handles.length_vessel - values_sec(n)).^2));
            ID_selected(n) = ID_nodes(I);
        end
        % [~,I] = min(handles.length_vessel(handles.length_vessel>0));
        min_value = min(handles.Laplace(handles.Laplace>0));
        I = ID_nodes(handles.Laplace==min_value);
        ID_selected(1) = ID_nodes(I(1));
        [~,I] = max(handles.length_vessel);
        ID_selected(end) = ID_nodes(I);
        
        % section of nodes.
        handles.NODES_SECTION = cell(length(ID_selected)-1,1);
        for n = 1:length(ID_selected)-1
            handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
        end
    
        % verifico que numero incluyo con un color y cual con otro color
        out_id = rem([1:length(ID_selected)-1], 2) == 0;
        
        % mostramos laplace en pantalla
        axes(handles.axes1);
        plot(0,0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
        hold on 
        for n = 1:length(out_id)
            if out_id(n)==0
                plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
            else
                plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
            end
        end
        plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
        plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
        hold off
        daspect([1,1,1])
        axis off
        view(handles.azimuth,handles.elevation)
    
    else
    
        
        ID_nodes = [1:size(handles.nodes,1)]';
        % laplace node selected
        ID_selected = zeros(handles.sec + 1,1);
    %     [~,I] = min(handles.length_vessel(handles.length_vessel>0));
    %     ID_selected(1) = ID_nodes(I);
        min_value = min(handles.Laplace(handles.Laplace>0));
        I = ID_nodes(handles.Laplace==min_value);
        ID_selected(1) = ID_nodes(I(1));
        
        [~,I] = max(handles.length_vessel);
        ID_selected(end) = ID_nodes(I);
        
        % section of nodes.
        handles.NODES_SECTION = cell(length(ID_selected)-1,1);
        for n = 1:length(ID_selected)-1
            handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
        end
    
        % verifico que numero incluyo con un color y cual con otro color
        out_id = rem([1:length(ID_selected)-1], 2) == 0;
        
        % mostramos laplace en pantalla
        axes(handles.axes1);
        plot(0,0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
        hold on 
        for n = 1:length(out_id)
            if out_id(n)==0
                plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
            else
                plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
            end
        end
        plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
        plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
        hold off
        daspect([1,1,1])
        axis off
        view(handles.azimuth,handles.elevation)
    end
end