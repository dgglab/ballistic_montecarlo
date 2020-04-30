function dev1 = runHallBar()

%device body properties
width = 2;
bodySQ = 2;
injSQ = 0.5;

    
%make the device body
dev1 = hallBar(width,bodySQ,injSQ);
dev1.plot()
pause(1)

%add ohmics
dev1.addEtch()
dev1.plot()
pause(1)


% pause(1)

%stitch overlaping features
dev1.joinFeatures();    
dev1.plot()
end