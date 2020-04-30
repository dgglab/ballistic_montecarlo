function dev1 = runDelafossiteBar(lengthBar,widthBar)

%device body properties
% lengthBar = 174;
% widthBar = 6.7;

    
%make the device body
dev1 = delafossiteBar(lengthBar,widthBar);
dev1.plot()
pause(1)

%add ohmics
dev1.addOhmics()
dev1.plot()
pause(1)


% pause(1)

%stitch overlaping features
dev1.joinFeatures();    
dev1.plot()
end