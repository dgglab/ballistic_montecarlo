function out = randcos(min_r,max_r)
    
            w_range=max_r-min_r;
            tempr=rand()*w_range+min_r;
            tempr1=rand();
            while tempr1>abs(cos(tempr))
               tempr =rand()*w_range+min_r;
               tempr1=rand();
            end
            
            out=tempr;