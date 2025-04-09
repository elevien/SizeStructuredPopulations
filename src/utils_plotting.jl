
function recursive_draw_tree(ax,node,x0,y0,depth,h,lw)
    val = 0. #max(0.2,min(abs(node.label/(τ*D)),1))
    clr = (val,val,val)
    dt = (node.label[end,1] -node.label[1,1])*h
    ax.plot([x0,x0],[y0,y0-dt],"-",c=clr,lw=lw)
    
    if node.left != nothing

        val =  0#max(0.2,min(abs(node.right.label/(τ*D)),1))
        clr = (val,val,val)
        ax.plot([x0,x0+2.0^(-depth)],[y0-dt,y0-dt],"-",c=clr,lw=lw)
        recursive_draw_tree(ax,node.left,x0-2.0^(-depth),y0-dt,depth+1,h,lw)

        val =  0#max(0.2,min(abs(node.left.label/(τ*D)),1))
        clr = (val,val,val)
        ax.plot([x0-2.0^(-depth),x0],[y0-dt,y0-dt],"-",c=clr,lw=lw)
        recursive_draw_tree(ax,node.right,x0+2.0^(-depth),y0-dt,depth+1,h,lw)
    end
end


function recursive_draw_lineage(ax,node,x0,y0,depth,h,lw,clr)

    dt = (node.label[end,1] -node.label[1,1])*h
    ax.plot([x0,x0],[y0,y0-dt],"-",c=clr,lw=lw)
    
    if node.left != nothing
        r = rand()
        if r <0.5
            ax.plot([x0,x0+2.0^(-depth)],[y0-dt,y0-dt],"-",c=clr,lw=lw)
            return recursive_draw_lineage(ax,node.right,x0+2.0^(-depth),y0-dt,depth+1,h,lw,clr)
        else
            ax.plot([x0-2.0^(-depth),x0],[y0-dt,y0-dt],"-",c=clr,lw=lw)
            return recursive_draw_lineage(ax,node.left,x0-2.0^(-depth),y0-dt,depth+1,h,lw,clr)
        end        
    end
    return x0,y0-dt,depth
end