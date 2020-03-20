classdef causticsBarrier < causticsFrame
    properties
         
        
       theta1;
       k_y;
       Tin
       Tout
       q_x
       thetain
       thetaout
       lastin=0
       E_rel
       
    end
    methods
        function obj = causticsBarrier(px,py,E_rel)
            obj = obj@causticsFrame(px,py,ones(1,length(px)-1));
            obj.E_rel=E_rel;
            obj.theta1=(-pi/2.00:.001:pi/2.00)';
            obj.k_y=sin(obj.theta1);
            r=obj.klein_tunneling_pn(obj.theta1,1,1-obj.E_rel);
            obj.Tin=real(1-r.*conj(r));

            obj.q_x=sign(obj.E_rel)*sqrt(obj.E_rel.^2-obj.k_y.^2);
            obj.thetain=real(atan(obj.k_y./obj.q_x));

            r=obj.klein_tunneling_pn(obj.theta1,1,-1./obj.E_rel+1);
            obj.Tout=real(1-r.*conj(r));
            obj.q_x=sign(1./obj.E_rel)*sqrt((1./obj.E_rel).^2-obj.k_y.^2);
            obj.thetaout=real(atan(obj.k_y./obj.q_x));

            
        end
        
        function theta_rel=rel_theta(obj,theta,crossed_line)
                       theta_rel=obj.phis(crossed_line)-pi/2-theta;
           if(cos(theta_rel)>0)
                theta_rel=mod(obj.phis(crossed_line)-theta,2*pi)-pi/2;
           else
                theta_rel=mod(obj.phis(crossed_line)+pi-theta,2*pi)-pi/2;
           end
           
        end
        
        function theta_rel=rel_theta2(obj,theta,crossed_line)
           theta_rel=mod(obj.edgenorm(crossed_line)-theta,2*pi);
           theta_rel(theta_rel>pi)=theta_rel(theta_rel>pi)-2*pi;
        end
        
        function out = reflect_off_barrier_in(obj,theta,crossed_line)
            theta_rel=rel_theta(obj,theta,crossed_line);
            out=(rand()>interp1(obj.theta1,obj.Tin,theta_rel));
        end                
               
        function theta_out = refract_in(obj,theta,crossed_line)
            theta_rel=rel_theta(obj,theta,crossed_line);
            theta_out=theta-interp1(obj.theta1,obj.thetain-obj.theta1,theta_rel);
        end
        
        function out = reflect_off_barrier_out(obj,theta,crossed_line)
            theta_rel=rel_theta(obj,theta,crossed_line);
            out=(rand()>interp1(obj.theta1,obj.Tout,theta_rel));
        end                
               
        function theta_out = refract_out(obj,theta,crossed_line)
            theta_rel=rel_theta(obj,theta,crossed_line);
            theta_out=theta-interp1(obj.theta1,obj.thetaout-obj.theta1,theta_rel);
        end
        
        function [theta_out,r_out]=reflect_or_transmit(obj,theta,r_global,crossed_line)
           if(abs(obj.rel_theta2(theta,crossed_line))<pi/2)
               if(obj.reflect_off_barrier_in(theta,crossed_line))
                   theta_out=obj.specular_reflection(theta,crossed_line);
                   r_out=r_global;
               else
                   theta_out=obj.refract_in(theta,crossed_line);
                   r_out=r_global*obj.E_rel;
               end
           else
               if(obj.reflect_off_barrier_out(theta,crossed_line))
                   theta_out=obj.specular_reflection(theta,crossed_line);
                   r_out=r_global*obj.E_rel;
               else
                   theta_out=obj.refract_out(theta,crossed_line);
                   r_out=r_global;
               end
           end
        end
        
        function [r,t]=klein_tunneling_pn(~,theta_in,E_f,dE_barrier)

            k_f=E_f;
            k_out=E_f-dE_barrier;

            k_x=k_f*cos(theta_in);
            k_y=k_f*sin(theta_in);


            s1=1;sign(k_f);
            s2=1;sign(k_out);
            q_x=s2.*sqrt(k_out.^2-k_y.^2);

            theta_out=atan(k_y./q_x);

            denom=1+exp(i*(theta_in+theta_out)).*s2.*s1;
            r=exp(i*theta_in).*(exp(i*theta_in)-exp(i*theta_out).*s1.*s2)./denom;
            t=(1+exp(2*i*theta_in))./denom;

        end
    end
end