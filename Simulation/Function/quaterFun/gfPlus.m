function result = gfPlus(x,y)
    gftable = [0,1,2,3;1,0,3,2;2,3,0,1;3,2,1,0];
    result = gftable(x+1,y+1);
end