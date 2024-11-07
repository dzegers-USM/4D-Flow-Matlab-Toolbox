function updateMainImageTexts(handles, c1, c2, b1, b2, a1, a2, inv)
    if inv == false
        set(handles.text1,'visible','on','String',['Slice: ',num2str(c1),' / ',num2str(c2)])
        set(handles.text116,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Type: ',handles.type])

        set(handles.text2,'visible','on','String',['Slice: ',num2str(b1),' / ',num2str(b2)])
        set(handles.text117,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Type: ',handles.type])

        set(handles.text3,'visible','on','String',['Slice: ',num2str(a1),' / ',num2str(a2)])
        set(handles.text118,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Type: ',handles.type])
    else
        set(handles.text1,'visible','on','String',['Slice: ',num2str(c1),' / ',num2str(c2)])
        set(handles.text116,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Type: ',handles.type])

        set(handles.text2,'visible','on','String',['Slice: ',num2str(b1),' / ',num2str(b2)])
        set(handles.text117,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Type: ',handles.type])

        set(handles.text3,'visible','on','String',['Slice: ',num2str(a1),' / ',num2str(a2)])
        set(handles.text118,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Type: ',handles.type])
    end
end
