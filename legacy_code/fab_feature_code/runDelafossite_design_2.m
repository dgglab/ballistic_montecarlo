function dev1 = runDelafossite_design_2

%device body properties
deviceLength = 21;
deviceWidth = 9;

% Top drains
dr1Width = 3;
dr1Pos = 2.5;
dr2Width = 3;
dr2Pos = 13.5;

% Injector array
inPos = 7.5; % Relative to left edge

    
%make the device body
dev1 = delafossite_design_2(deviceLength,deviceWidth);
dev1.plot()
pause(1)

%add drains
dev1.addDrains(dr1Width,dr1Pos,dr2Width,dr2Pos)
dev1.plot()
pause(1)

%add injectors
dev1.addInjectors(inPos)
dev1.plot()
pause(1)

% pause(1)

%stitch overlaping features
dev1.joinFeatures();    
dev1.plot()
end