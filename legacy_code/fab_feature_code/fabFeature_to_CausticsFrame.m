function [causticsFrameOut,ohmicsframes,barrierframes]= fabFeature_to_CausticsFrame(fabFeatureIn)

    ohmics_idx=[];
    barrier_idx=[];
    ohmicsframes={};
    barrierframes={};
    for i=1:fabFeatureIn.N_Features
        if fabFeatureIn.fabFeatures{i}.LayerNum==0
            frame.px=fabFeatureIn.fabFeatures{i}.px;
            frame.py=fabFeatureIn.fabFeatures{i}.py;
        elseif fabFeatureIn.fabFeatures{i}.LayerNum==1
            ohmics_idx=[ohmics_idx, i];
        elseif fabFeatureIn.fabFeatures{i}.LayerNum==2
            barrier_idx=[barrier_idx, i];
        end
        
    end
    
    for i=1:length(ohmics_idx)
        ohmicsframes{i}.px=fabFeatureIn.fabFeatures{ohmics_idx(i)}.px;
        ohmicsframes{i}.py=fabFeatureIn.fabFeatures{ohmics_idx(i)}.py;
    end
    for i=1:length(barrier_idx)
        barrierframes{i}.px=fabFeatureIn.fabFeatures{barrier_idx(i)}.px;
        barrierframes{i}.py=fabFeatureIn.fabFeatures{barrier_idx(i)}.py;
    end

    causticsFrameOut=causticsFrame(frame.px,frame.py,zeros(size(frame.px)));
    causticsFrameOut.add_ohmics(ohmicsframes);

end