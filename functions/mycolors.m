function col = mycolors(var,n)
    switch var
        case 'Mf'
            col = [1,0,0;...
                   0,.6,0;...
                   0,0,1;...
                   .7,0,.7];
        case 'Md'
            col = lines(7);
            col = col([3,2,5,1,7],:);
    end
end