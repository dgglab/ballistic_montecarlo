function plotFrameGroup(frmgrp_in)

    set(gcf,'color','white');
    set(gca,'visible','off');
    Nfrms=length(frmgrp_in.frms);
    for k=1:Nfrms
          switch class(frmgrp_in.frms{k})
               case 'causticsFrame'
                   frm1=frmgrp_in.frms{k};
                    plot(frm1.px,frm1.py,'k','linewidth',1.5);
                    hold on;
                    for i=1:length(frm1.ohmics)
                        if sum(frm1.terminalOhmics==frm1.ohmics{i}.ohmNumber)
                            plot(frm1.ohmics{i}.px,frm1.ohmics{i}.py,'k','Linewidth',2);
                        elseif sum(frm1.injector_ohmic==frm1.ohmics{i}.ohmNumber)
                            plot(frm1.ohmics{i}.px,frm1.ohmics{i}.py,'r','Linewidth',2);
                        else
                            plot(frm1.ohmics{i}.px,frm1.ohmics{i}.py,'b','Linewidth',2);
                        end
                    end
              case 'causticsBarrier'
                  cb1=frmgrp_in.frms{k};
                  plot(cb1.px,cb1.py,'color',[0,.5,0]);
                  hold on;
          
          end
          
    
    end
    hold off
    axis equal
    
    
    