function icon_label = labelWithIcon(label, icon_loc)
    if isdeployed
        root = strcat([ctfroot, '/APP_4D_FLOW']);
    else
        root = pwd;
    end
    root = replace(root, "\", "/");
    file_loc = strcat(['file:/', root, icon_loc]);
    icon_label = strcat(['<html><img width=16 height=16 src="', file_loc, '"></img>', label]);
end
