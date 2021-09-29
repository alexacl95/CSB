function res = unity(Uint, Lint, range)

    Uint = filloutliers(Uint,'center',2);
    Lint = filloutliers(Lint,'center',2);

    Uint = (Uint-range(:,1))./(range(:,2)-range(:,1));
    Lint = (Lint-range(:,1))./(range(:,2)-range(:,1));


    Uint = sort(Uint,2,'ascend');
    Lint = sort(Lint,2,'ascend');

    [Np,~] = size(Uint);
    
%     len = length(Uint);
%     Uint = Uint(:,round(0.25*len)+1:round(0.75*len));
%     Lint = Lint(:,round(0.25*len)+1:round(0.75*len));
    
    maxU = max(Uint,[],2);
    minU = min(Uint,[],2);
    maxL = max(Lint,[],2);
    minL = min(Lint,[],2);

    res = 1-(sum(maxU-minU) + sum(maxL-minL))/(Np*2);

end