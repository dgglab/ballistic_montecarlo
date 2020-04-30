

norms = deln.unorms;
norms = [0,-pi];
indstore = zeros(1E6,length(norms));
for edgeind = 1:length(norms)
    
    
    for i = 1:length(indstore)
        indstore(i,edgeind) = deln.inject(norms(edgeind));
    end
    
    figure(edgeind); clf;
    histogram(indstore(:,edgeind),1:253)
    xlim([0,253])
end

figure(6)
histogram(squeeze(indstore),0:255)
xlim([0,253])