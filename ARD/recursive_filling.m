function [M] = recursive_filling(vec,M,dir,idcol,idrow,maxcoll,maxrow,minx,miny)

    if(size(vec,2) == 0)
        return;
    end

    if(strcmp(dir,'down'))
        if (idrow<maxrow)
            M(idrow,idcol) =  vec(1);
            M = recursive_filling(vec(1,2:end),M,dir,idcol,idrow+1,maxcoll,maxrow,minx,miny);
        else
            M = recursive_filling(vec,M,'left',idcol-1,idrow-1,maxcoll-1,maxrow,minx,miny);
        end
    elseif(strcmp(dir,'left'))
        
        if (idcol>minx)
            M(idrow,idcol) =  vec(1);
            M = recursive_filling(vec(1,2:end),M,dir,idcol-1,idrow,maxcoll,maxrow,minx,miny);
        else
            M = recursive_filling(vec,M,'up',idcol+1,idrow-1,maxcoll,maxrow-1,minx,miny);
        end
    elseif(strcmp(dir,'up'))
        if (idrow>miny)
            M(idrow,idcol) =  vec(1);
            M = recursive_filling(vec(1,2:end),M,dir,idcol,idrow-1,maxcoll,maxrow,minx,miny);
        else
            M = recursive_filling(vec,M,'right',idcol+1,idrow+1,maxcoll,maxrow,minx+1,miny);
        end
    else %right
        if (idcol<maxcoll)
            M(idrow,idcol) =  vec(1);
            M = recursive_filling(vec(1,2:end),M,dir,idcol+1,idrow,maxcoll,maxrow,minx,miny);
        else
            M = recursive_filling(vec,M,'down',idcol-1,idrow+1,maxcoll,maxrow,minx,miny+1);
        end
    end
end